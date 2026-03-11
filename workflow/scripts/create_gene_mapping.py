#!/usr/bin/env python3

import os
import re
import argparse
import sys

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Create gene ID mapping table from merged GTF file'
    )
    parser.add_argument('--merged-gtf', required=True,
                        help='Merged GTF file from StringTie')
    parser.add_argument('--output', required=True,
                        help='Output mapping file (TSV format)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    return parser.parse_args()

def extract_gene_mapping(merged_gtf_file, verbose=False):
    """Extract gene ID mapping from merged GTF file"""
    gene_mapping = {}
    
    with open(merged_gtf_file) as f:
        for line_no, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            # Only process transcript lines
            if fields[2] != 'transcript':
                continue
            
            attributes = fields[8]
            
            # Extract gene_id and ref_gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            ref_gene_id_match = re.search(r'ref_gene_id "([^"]+)"', attributes)
            
            if gene_id_match and ref_gene_id_match:
                stringtie_id = gene_id_match.group(1)
                original_id = ref_gene_id_match.group(1)
                
                if stringtie_id in gene_mapping:
                    if gene_mapping[stringtie_id] != original_id:
                        print(f"Warning: Conflicting mapping for {stringtie_id}: "
                              f"{gene_mapping[stringtie_id]} vs {original_id}",
                              file=sys.stderr)
                else:
                    gene_mapping[stringtie_id] = original_id
                    
                    if verbose and len(gene_mapping) <= 10:
                        print(f"Mapping: {stringtie_id} -> {original_id}")
    
    return gene_mapping

def main():
    args = parse_arguments()

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    if args.verbose:
        print(f"Reading merged GTF file: {args.merged_gtf}")
    
    # Extract gene mapping
    gene_mapping = extract_gene_mapping(args.merged_gtf, args.verbose)
    
    if args.verbose:
        print(f"Found {len(gene_mapping)} gene mappings")
    
    # Write mapping to file
    with open(args.output, 'w') as f:
        f.write("StringTie_ID\tOriginal_ID\n")
        for stringtie_id, original_id in sorted(gene_mapping.items()):
            f.write(f"{stringtie_id}\t{original_id}\n")
    
    print(f"Gene mapping written to: {args.output}")
    print(f"Total mappings: {len(gene_mapping)}")

if __name__ == "__main__":
    main()
