# Bombyx mori Non-model Insect Case Study — Pipeline Configuration

## Dataset
- **Source**: DRA008737 (P50T strain, 5th instar larvae)
- **Tissues**: Testis (TT, male) vs Ovary (OV, female)
- **Contrast**: 3 male vs 3 female
- **SRA Runs**: DRR186498-DRR186500 (TT), DRR186501-DRR186503 (OV)

## Download Commands
```bash
conda activate kingfisher

# Testis (male) — 3 replicates
kingfisher get -r DRR186498 -m ena-ftp -o data/fastq/bombyx/
kingfisher get -r DRR186499 -m ena-ftp -o data/fastq/bombyx/
kingfisher get -r DRR186500 -m ena-ftp -o data/fastq/bombyx/

# Ovary (female) — 3 replicates
kingfisher get -r DRR186501 -m ena-ftp -o data/fastq/bombyx/
kingfisher get -r DRR186502 -m ena-ftp -o data/fastq/bombyx/
kingfisher get -r DRR186503 -m ena-ftp -o data/fastq/bombyx/
```

## Reference Genome
- **Assembly**: Kawamoto et al. 2019
- **Annotation**: GFF3 format (available alongside the genome)
- Place genome FASTA in `data/reference/bombyx/genome.fasta`
- Place annotation in `data/reference/bombyx/annotation.gff3`

## After Download — Run Pipeline
```bash
# 1. Update config.yaml:
#    reference.genome: "data/reference/bombyx/genome.fasta"
#    reference.annotation: "data/reference/bombyx/annotation.gff3"
#    samples: "data/fastq/bombyx/samples.tsv"

# 2. Create data/fastq/bombyx/samples.tsv:
#    sample	fq1	fq2	group
#    TT_1	DRR186498_1.fastq.gz	DRR186498_2.fastq.gz	Testis
#    TT_2	DRR186499_1.fastq.gz	DRR186499_2.fastq.gz	Testis
#    TT_3	DRR186500_1.fastq.gz	DRR186500_2.fastq.gz	Testis
#    OV_1	DRR186501_1.fastq.gz	DRR186501_2.fastq.gz	Ovary
#    OV_2	DRR186502_1.fastq.gz	DRR186502_2.fastq.gz	Ovary
#    OV_3	DRR186503_1.fastq.gz	DRR186503_2.fastq.gz	Ovary

# 3. dea.comparisons: "Testis_vs_Ovary"

# 4. Run pipeline:
#    ./run_analysis.sh
```

## Expected Case Study Outputs
After pipeline completion, verify:
1. Reference auto-detection correctly identified non-standard file names
2. Namespace module resolved Bombyx-specific gene IDs  
3. Consensus produced Tier A/B/C genes
4. Decontam audit shows real-world contamination levels (if any)
5. Functional enrichment of Tier A genes is biologically plausible
