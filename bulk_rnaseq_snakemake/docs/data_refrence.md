# Reference and Raw Data Preparation

## Project

GSE334005 Mouse Bulk RNA-seq

- Organism: Mus musculus
- Genome assembly: GRCm39
- RNA-seq layout: Paired-end
- Number of samples: 30
- Aligner: HISAT2
- Annotation: GENCODE M39 / Ensembl 116
- Workflow manager: Snakemake

---

# 1. Create Project Directories

Go to the project directory:

```bash
cd ~/bulk_rnaseq_snakemake


bulk_rnaseq_snakemake/
│
├── data/
│   ├── fastq/
│   └── sra/
│
├── reference/
│   ├── GRCm39.primary_assembly.genome.fa
│   ├── gencode.vM39.primary_assembly.annotation.gtf
│   └── hisat2_index/
│
├── metadata/
│
├── scripts/
│
└── results/
    ├── fastqc/
    ├── trimmed/
    ├── aligned/
    ├── counts/
    ├── deseq2/
    ├── GO/
    └── KEGG/



