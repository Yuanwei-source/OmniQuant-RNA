import os
import sys
import glob
import re
import argparse

def find_fastq_files(directory):
    """Find all FASTQ files in the directory."""
    fastq_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    files = []
    for root, _, filenames in os.walk(directory):
        for filename in filenames:
            if any(filename.endswith(ext) for ext in fastq_extensions):
                files.append(os.path.join(root, filename))
    return files

def extract_sample_info(files):
    """Extract sample names and paths from filenames."""
    samples = {}
    for file_path in files:
        filename = os.path.basename(file_path)
        # Remove extension
        name = filename
        for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        
        # Remove read suffix (_R1, _1, etc.)
        name = re.sub(r'[_.-]?(R?)[12]$', '', name)
        
        # Store absolute path. Prefer R1 if available, or just the first one found.
        # If we already have this sample, check if the new file is R1 (usually better for "path")
        abs_path = os.path.abspath(file_path)
        
        if name not in samples:
            samples[name] = abs_path
        else:
            # If current stored path is R2 and new is R1, swap it?
            # Simple heuristic: if "R1" or "_1" in new path and not in old, update.
            if ("R1" in abs_path or "_1." in abs_path) and ("R1" not in samples[name] and "_1." not in samples[name]):
                samples[name] = abs_path
                
    return samples

def main():
    parser = argparse.ArgumentParser(
        description='Generate sample configuration from FASTQ files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate sample list
  python generate_samples.py data/fastq/
  
  # Output to file
  python generate_samples.py data/fastq/ -o samples.tsv
        """
    )
    
    parser.add_argument('fastq_dir', help='Directory containing FASTQ files')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.fastq_dir):
        print(f"Error: {args.fastq_dir} is not a directory")
        sys.exit(1)
    
    # Find FASTQ files
    fastq_files = find_fastq_files(args.fastq_dir)
    if not fastq_files:
        print(f"No FASTQ files found in {args.fastq_dir}")
        sys.exit(1)
    
    # Extract sample info
    samples_dict = extract_sample_info(fastq_files)
    
    # Prepare output
    output_lines = []
    
    # Header
    header = ["sample", "path"]
    output_lines.append('\t'.join(header))
    
    # Data rows
    for sample in sorted(samples_dict.keys()):
        output_lines.append(f"{sample}\t{samples_dict[sample]}")
    
    # Output
    if args.output:
        with open(args.output, 'w') as f:
            f.write('\n'.join(output_lines) + '\n')
        print(f"Sample configuration written to {args.output}")
    else:
        for line in output_lines:
            print(line)

if __name__ == "__main__":
    main()
