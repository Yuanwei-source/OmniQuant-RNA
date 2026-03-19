#!/usr/bin/env python3

import os
import sys
import argparse
import polars as pl
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='Aggregate StringTie results using native ctab/tab files')
    parser.add_argument('-i', '--input_dir', required=True, help="Base directory of sample results")
    parser.add_argument('-p', '--pattern', help="Pattern for sample files (ignored, using ctab automatically)")
    parser.add_argument('--gene-mapping', help="Gene ID mapping file")
    parser.add_argument('--output-gene-counts', required=True)
    parser.add_argument('--output-transcript-counts', required=True)
    parser.add_argument('--output-gene-tpm', required=True)
    parser.add_argument('--output-transcript-tpm', required=True)
    parser.add_argument('--output-gene-fpkm', help="Optional gene FPKM output")
    parser.add_argument('--output-transcript-fpkm', help="Optional transcript FPKM output")
    parser.add_argument('-l', '--length', type=int, default=150, help="Read length for count estimation")
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser.parse_args()

def load_mapping(path):
    if not path or not os.path.exists(path): return {}
    df = pl.read_csv(path, separator='\t')
    if len(df.columns) < 2: return {}
    return dict(zip(df.to_series(0), df.to_series(1)))

def process_samples(args):
    samples = []
    # Find all t_data.ctab within the input_dir
    for ctab_path in Path(args.input_dir).rglob("t_data.ctab"):
        parts = ctab_path.parts
        if "final" in parts:
            sample_name = parts[parts.index("final") - 1]
            samples.append((sample_name, ctab_path.parent))

    gene_mapping = load_mapping(args.gene_mapping)
    
    t_counts_list, t_tpm_list, t_fpkm_list = [], [], []
    g_tpm_list, g_fpkm_list = [], []
    sample_names = []

    for name, s_dir in sorted(samples):
        sample_names.append(name)
        
        # 1. Handle Transcripts (t_data.ctab)
        ctab = pl.read_csv(s_dir / "t_data.ctab", separator='\t')
        
        # coverage to counts formula as done by stringtie's prepDE
        ctab = ctab.with_columns([
            ((pl.col("cov") * pl.col("length")) / args.length).ceil().cast(pl.Int64).alias("counts"),
            (pl.col("cov") * 10).alias("TPM") # TPM in stringtie t_data.ctab is approximated by cov? Actually StringTie does not output TPM in t_data.ctab. Wait! Stringtie t_data.ctab format: t_id, chr, strand, start, end, t_name, num_exons, length, gene_id, gene_name, cov, FPKM. Let's just output cov and FPKM if TPM is unavailable, or calculate TPM from FPKM. TPM = (FPKM / sum(FPKM)) * 1e6
        ])
        
        # calculate TPM accurately if needed:
        sum_fpkm = ctab["FPKM"].sum()
        if sum_fpkm > 0:
            ctab = ctab.with_columns([
                ((pl.col("FPKM") / sum_fpkm) * 1e6).alias("TPM_calc")
            ])
        else:
            ctab = ctab.with_columns(pl.lit(0.0).alias("TPM_calc"))
            
        t_counts_list.append(ctab.select(["t_name", "counts"]).rename({"counts": name}))
        t_tpm_list.append(ctab.select(["t_name", "TPM_calc"]).rename({"TPM_calc": name}))
        t_fpkm_list.append(ctab.select(["t_name", "FPKM"]).rename({"FPKM": name}))

        # 2. Handle Genes (gene_abundances.tab)
        gtab = pl.read_csv(s_dir / "gene_abundances.tab", separator='\t', comment_prefix='#', ignore_errors=True)
        gid_col = "Gene ID" if "Gene ID" in gtab.columns else gtab.columns[0]
        
        tpm_col = "TPM" if "TPM" in gtab.columns else None
        fpkm_col = "FPKM" if "FPKM" in gtab.columns else None
        
        if gene_mapping:
            gtab = gtab.with_columns(
                pl.col(gid_col).replace(gene_mapping, default=pl.col(gid_col)).alias("Gene_ID")
            )
        else:
            gtab = gtab.rename({gid_col: "Gene_ID"})
            
        if tpm_col:
            g_tpm_list.append(gtab.select(["Gene_ID", tpm_col]).rename({tpm_col: name}))
        if fpkm_col:
            g_fpkm_list.append(gtab.select(["Gene_ID", fpkm_col]).rename({fpkm_col: name}))

    if not samples:
        print("No t_data.ctab found.", file=sys.stderr)
        return

    def join_all(df_list, join_key):
        base = df_list[0]
        for next_df in df_list[1:]:
            try:
                base = base.join(next_df, on=join_key, how="outer_coalesce").fill_null(0)
            except Exception:
                base = base.join(next_df, on=join_key, how="outer").fill_null(0)
                # coalesce manually if duplicated keys exist
                if join_key + "_right" in base.columns:
                    base = base.with_columns(
                        pl.col(join_key).fill_null(pl.col(join_key + "_right"))
                    ).drop(join_key + "_right")
        return base

    t_counts_joined = join_all(t_counts_list, "t_name")
    t_counts_joined.write_csv(args.output_transcript_counts, separator='\t')
    
    t_tpm_joined = join_all(t_tpm_list, "t_name")
    t_tpm_joined.write_csv(args.output_transcript_tpm, separator='\t')
    
    if args.output_transcript_fpkm:
        t_fpkm_joined = join_all(t_fpkm_list, "t_name")
        t_fpkm_joined.write_csv(args.output_transcript_fpkm, separator='\t')

    if g_tpm_list:
        join_all(g_tpm_list, "Gene_ID").write_csv(args.output_gene_tpm, separator='\t')
    if args.output_gene_fpkm and g_fpkm_list:
        join_all(g_fpkm_list, "Gene_ID").write_csv(args.output_gene_fpkm, separator='\t')

    # Gene Counts Aggregation using transcript mapping
    # t_data.ctab has 'gene_id'
    t2g = ctab.select(["t_name", "gene_id"])
    if gene_mapping:
        t2g = t2g.with_columns(pl.col("gene_id").replace(gene_mapping, default=pl.col("gene_id")))
    
    # join and group by gene_id
    gene_counts = t2g.join(t_counts_joined, on="t_name").drop("t_name")
    
    # Sum counts per gene_id
    gene_counts = gene_counts.group_by("gene_id").sum()
    gene_counts.rename({"gene_id": "Gene_ID"}).write_csv(args.output_gene_counts, separator='\t')

if __name__ == "__main__":
    process_samples(parse_arguments())
