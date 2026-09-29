# Bulk RNA-seq Workflow

## Project

This repository contains an automated bulk RNA-seq analysis pipeline for mouse RNA-seq data from GEO accession:

**GSE334005**

The dataset contains:

- Organism: *Mus musculus*
- Number of samples: 30
- Sequencing type: Paired-end
- Reference genome: GRCm39
- Annotation: GENCODE M39 / Ensembl 116
- Aligner: HISAT2
- Quantification: featureCounts
- Differential expression: DESeq2
- Functional enrichment: GO and KEGG
- Workflow manager: Snakemake

---

# 1. Overall Workflow

The complete analysis follows:

```text
GEO / SRA
   │
   ▼
raw_data.sh
   │
   ├── Download 30 SRA samples
   ├── Convert SRA → paired FASTQ
   └── Prepare reference genome
            │
            └── HISAT2 index
   │
   ▼
Snakemake
   │
   ├── FastQC
   │
   ├── MultiQC
   │
   ├── fastp
   │
   ├── HISAT2
   │
   ├── SAM/BAM processing
   │
   ├── featureCounts
   │
   ├── Count matrix
   │
   ├── DESeq2
   │
   ├── GO enrichment
   │
   └── KEGG enrichment
   │
   ▼
Final Results
