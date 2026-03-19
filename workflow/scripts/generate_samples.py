import os
import sys
import glob
import re
import argparse
from pathlib import Path

READ1_PATTERNS = [re.compile(p, flags=re.IGNORECASE) for p in [r'_R1(?:\D|$)', r'\.R1(?:\D|$)', r'_1(?:\D|$)', r'\.1(?:\D|$)']]
READ2_PATTERNS = [re.compile(p, flags=re.IGNORECASE) for p in [r'_R2(?:\D|$)', r'\.R2(?:\D|$)', r'_2(?:\D|$)', r'\.2(?:\D|$)']]

def find_fastq_files(directory):
    """Find all FASTQ files in the directory."""
    fastq_extensions = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}
    dir_path = Path(directory)
    files = []
    for p in dir_path.rglob('*'):
        if p.is_file() and any(p.name.endswith(ext) for ext in fastq_extensions):
            files.append(str(p))
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
        # Handle .R1.fastq.gz pattern
        base_name = name
        if '.R1' in base_name:
            base_name = base_name.replace('.R1', '')
        elif '.R2' in base_name:
            base_name = base_name.replace('.R2', '')
        else:
            base_name = re.sub(r'[_.-]?(R?)[12]$', '', base_name)
        
        if base_name not in samples:
            samples[base_name] = {'r1': None, 'r2': None}
            
        abs_path = os.path.abspath(file_path)
        
        filename_lower = filename.lower()
        is_r1 = any(p.search(filename_lower) for p in READ1_PATTERNS)
        is_r2 = any(p.search(filename_lower) for p in READ2_PATTERNS)

        if is_r1 and is_r2:
            raise ValueError(f"Ambiguous read direction for file: {file_path}")

        if is_r1:
            if samples[base_name]['r1'] is not None:
                raise ValueError(f"Duplicate R1 for sample {base_name}: {samples[base_name]['r1']} and {abs_path}")
            samples[base_name]['r1'] = abs_path
        elif is_r2:
            if samples[base_name]['r2'] is not None:
                raise ValueError(f"Duplicate R2 for sample {base_name}: {samples[base_name]['r2']} and {abs_path}")
            samples[base_name]['r2'] = abs_path
        else:
            raise ValueError(
                f"Cannot infer read direction (R1/R2) from file name: {file_path}. "
                "Please use *_R1/*_R2 or *_1/*_2 naming."
            )
                
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
    try:
        samples_dict = extract_sample_info(fastq_files)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    
    # Prepare output
    output_lines = []
    
    # Header
    header = ["sample", "fq1", "fq2", "group"]
    output_lines.append('\t'.join(header))
    
    # Data rows
    for sample in sorted(samples_dict.keys()):
        r1 = samples_dict[sample]['r1']
        r2 = samples_dict[sample]['r2']
        if r1:
            line = f"{sample}\t{r1}\t{r2 if r2 else ''}\t"
            output_lines.append(line)
    
    # Output
    if args.output:
        output_path = args.output
        # If output path is an existing directory, append samples.tsv
        if os.path.isdir(output_path):
            output_path = os.path.join(output_path, "samples.tsv")
            
        with open(output_path, 'w') as f:
            f.write('\n'.join(output_lines) + '\n')
        print(f"Sample configuration written to {output_path}")
    else:
        for line in output_lines:
            print(line)

if __name__ == "__main__":
    main()
