import re

filenames = [
    "Snakefile",
    "workflow/rules/quantification_stringtie.smk",
    ".gitignore"
]

def update_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    # We want to change the raw string "gene_counts_matrix.tsv" to "results/04.quantification/stringtie/gene_counts_matrix.tsv" etc.
    replacements = {
        '"gene_counts_matrix.tsv"': '"results/04.quantification/stringtie/gene_counts_matrix.tsv"',
        '"transcript_counts_matrix.tsv"': '"results/04.quantification/stringtie/transcript_counts_matrix.tsv"',
        '"gene_tpm_matrix.tsv"': '"results/04.quantification/stringtie/gene_tpm_matrix.tsv"',
        '"transcript_tpm_matrix.tsv"': '"results/04.quantification/stringtie/transcript_tpm_matrix.tsv"',
    }
    
    # The only exception is .gitignore, let's just make sure it stays ignored without prepending quotes.
    if "gitignore" in filepath:
        content = content.replace("gene_counts_matrix.tsv", "results/04.quantification/stringtie/gene_counts_matrix.tsv")
        content = content.replace("gene_tpm_matrix.tsv", "results/04.quantification/stringtie/gene_tpm_matrix.tsv")
        content = content.replace("transcript_counts_matrix.tsv", "results/04.quantification/stringtie/transcript_counts_matrix.tsv")
        content = content.replace("transcript_tpm_matrix.tsv", "results/04.quantification/stringtie/transcript_tpm_matrix.tsv")
    else:
        for old, new in replacements.items():
            content = content.replace(old, new)
            
    with open(filepath, 'w') as f:
        f.write(content)

for fname in filenames:
    update_file(fname)

print("Done")
