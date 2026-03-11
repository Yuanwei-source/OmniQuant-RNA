#!/usr/bin/env python3
"""Build sample-to-file import manifests for tximport-based quantifiers."""

import argparse
import csv
from pathlib import Path

import pandas as pd


def build_manifest(samples, quantifier, import_type, file_role, resolver):
    rows = []
    for sample in samples:
        file_path = resolver(sample)
        if not Path(file_path).exists():
            raise FileNotFoundError(f"Missing expected {quantifier} file for sample {sample}: {file_path}")
        rows.append(
            {
                "sample": sample,
                "file_path": file_path,
                "quantifier": quantifier,
                "import_type": import_type,
                "file_role": file_role,
            }
        )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Build tximport manifests")
    parser.add_argument("--samples", required=True, help="Sample sheet TSV")
    parser.add_argument("--salmon-output", required=True, help="Output path for salmon manifest")
    parser.add_argument("--kallisto-output", required=True, help="Output path for kallisto manifest")
    parser.add_argument("--stringtie-output", required=True, help="Output path for stringtie manifest")
    args = parser.parse_args()

    samples_df = pd.read_csv(args.samples, sep="\t")
    samples = samples_df["sample"].tolist()

    salmon = build_manifest(
        samples,
        quantifier="salmon",
        import_type="salmon",
        file_role="quant.sf",
        resolver=lambda sample: f"results/04.quantification/salmon/{sample}/quant.sf",
    )
    kallisto = build_manifest(
        samples,
        quantifier="kallisto",
        import_type="kallisto",
        file_role="abundance.tsv",
        resolver=lambda sample: f"results/04.quantification/kallisto/{sample}/abundance.tsv",
    )
    stringtie = build_manifest(
        samples,
        quantifier="stringtie",
        import_type="stringtie",
        file_role="t_data.ctab",
        resolver=lambda sample: f"results/04.quantification/stringtie/{sample}/final/t_data.ctab",
    )

    for path, df in [
        (args.salmon_output, salmon),
        (args.kallisto_output, kallisto),
        (args.stringtie_output, stringtie),
    ]:
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
        print(f"[import_manifest] Wrote {len(df)} rows to {output_path}")


if __name__ == "__main__":
    main()
