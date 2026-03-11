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
        
        # Determine if R1 or R2
        if '_R1' in filename or '_1.' in filename or filename.endswith('_1') or '.R1.' in filename:
            samples[base_name]['r1'] = abs_path
        elif '_R2' in filename or '_2.' in filename or filename.endswith('_2') or '.R2.' in filename:
            samples[base_name]['r2'] = abs_path
        else:
            # Fallback: if we can't determine, assume R1 if empty
            if samples[base_name]['r1'] is None:
                samples[base_name]['r1'] = abs_path
                
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
        print("Error: {} is not a directory".format(args.fastq_dir))
        sys.exit(1)
    
    # Find FASTQ files
    fastq_files = find_fastq_files(args.fastq_dir)
    if not fastq_files:
        print("No FASTQ files found in {}".format(args.fastq_dir))
        sys.exit(1)
    
    # Extract sample info
    samples_dict = extract_sample_info(fastq_files)
    
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
            line = "{}\t{}\t{}\t".format(sample, r1, r2 if r2 else '')
            output_lines.append(line)
    
    # Output
    if args.output:
        output_path = args.output
        # If output path is an existing directory, append samples.tsv
        if os.path.isdir(output_path):
            output_path = os.path.join(output_path, "samples.tsv")
            
        with open(output_path, 'w') as f:
            f.write('\n'.join(output_lines) + '\n')
        print("Sample configuration written to {}".format(output_path))
    else:
        for line in output_lines:
            print(line)

if __name__ == "__main__":
    main()
