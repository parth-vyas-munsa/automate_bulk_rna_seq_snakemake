# Automated Mouse Bulk RNA-seq Analysis Pipeline

## 📌 Project Overview

This repository contains a **fully automated, reproducible bulk RNA-seq analysis pipeline** for mouse RNA-seq data using **Snakemake**.

The pipeline starts from publicly available **SRA sequencing data** and processes the samples through quality control, adapter/quality trimming, genome alignment, BAM processing, read quantification, differential expression analysis, and functional enrichment.

The project was developed as a practical bioinformatics workflow demonstrating how a multi-sample RNA-seq analysis can be automated using a workflow management system rather than executing every command manually.

---

## 🧬 Dataset

**GEO accession:** `GSE334005`

**Organism:** *Mus musculus*

**Genome:** GRCm39

**Sequencing:** Paired-end RNA-seq

**Total samples:** 30

### Experimental groups

The dataset contains three treatment groups:

| Treatment | Number of samples |
| --------- | ----------------: |
| PBS       |                10 |
| OVM       |                10 |
| OVM18     |                10 |
| **Total** |            **30** |

The samples also contain three cell lines:

| Cell line |    PBS |    OVM |  OVM18 |  Total |
| --------- | -----: | -----: | -----: | -----: |
| CT26      |      4 |      4 |      4 |     12 |
| B16       |      3 |      3 |      3 |      9 |
| LLC       |      3 |      3 |      3 |      9 |
| **Total** | **10** | **10** | **10** | **30** |

The metadata are stored separately and are used throughout the downstream analysis.

---

# 🔬 Workflow

The complete workflow is:

```text
SRA / GEO
    │
    ▼
Download SRA accessions
    │
    ▼
Convert SRA → paired FASTQ
    │
    ▼
FastQC
    │
    ▼
MultiQC
    │
    ▼
fastp
    │
    ▼
HISAT2 alignment
    │
    ▼
SAM → BAM
    │
    ▼
Sort BAM
    │
    ▼
Index BAM
    │
    ▼
featureCounts
    │
    ▼
Count Matrix
    │
    ▼
DESeq2
    │
    ├── OVM vs PBS
    ├── OVM18 vs PBS
    └── OVM18 vs OVM
    │
    ▼
Significant DEGs
    │
    ├── Upregulated genes
    └── Downregulated genes
    │
    ▼
GO enrichment
    │
    ▼
KEGG enrichment
    │
    ▼
Figures + Tables + Results
```

---

# 🛠️ Tools Used

| Tool            | Purpose                            |
| --------------- | ---------------------------------- |
| SRA Toolkit     | Download and convert SRA data      |
| FastQC          | Raw read quality control           |
| MultiQC         | Aggregate QC reports               |
| fastp           | Adapter and quality trimming       |
| HISAT2          | Read alignment                     |
| SAMtools        | BAM processing                     |
| featureCounts   | Gene-level read quantification     |
| Snakemake       | Workflow automation                |
| R               | Statistical analysis               |
| DESeq2          | Differential expression            |
| clusterProfiler | GO/KEGG enrichment                 |
| org.Mm.eg.db    | Mouse gene annotation              |
| ggplot2         | Visualization                      |
| pheatmap        | Heatmaps                           |
| ggrepel         | Gene labeling                      |
| Python/Bash     | Workflow and preprocessing scripts |

---

# 📂 Repository Structure

```text
bulk-rnaseq-snakemake/
│
├── README.md
├── Snakefile
├── config.yaml
├── environment.yml
├── .gitignore
│
├── metadata/
│   └── metadata.csv
│
├── scripts/
│   ├── raw_data.sh
│   ├── merge_counts.R
│   ├── deseq2.R
│   └── go_kegg.R
│
├── docs/
│   ├── REFERENCE.md
│   ├── WORKFLOW.md
│   └── DESEQ2_CONNECT.md
│
├── reference/
│   ├── Mus_musculus.GRCm39.fa
│   ├── Mus_musculus.GRCm39.116.gtf
│   └── hisat2_index/
│
├── data/
│   └── fastq/
│
└── results/
    ├── fastqc/
    ├── multiqc/
    ├── trimmed/
    ├── aligned/
    ├── counts/
    ├── deseq2/
    └── enrichment/
```

> Large sequencing and reference files are intentionally excluded from GitHub using `.gitignore`.

---

# 1. Download Raw Data

The script:

```text
scripts/raw_data.sh
```

handles the initial project setup.

It performs tasks such as:

* Creating required directories
* Downloading SRA accessions
* Converting SRA files to FASTQ
* Downloading the mouse reference genome
* Downloading the annotation file
* Building the HISAT2 reference index

Example:

```bash
bash scripts/raw_data.sh
```

After execution, the project should contain the required input data and reference files.

---

# 2. Quality Control

## FastQC

FastQC is used to evaluate the quality of raw sequencing reads.

Typical metrics include:

* Per-base sequence quality
* Per-sequence quality
* Adapter content
* Sequence duplication
* GC content
* Overrepresented sequences
* Sequence length

The results are stored in:

```text
results/fastqc/
```

---

# 3. MultiQC

MultiQC combines the individual FastQC reports into one summary report.

Output:

```text
results/multiqc/multiqc_report.html
```

This makes it easier to inspect all 30 samples together.

---

# 4. Read Trimming

The pipeline uses **fastp** for preprocessing.

fastp performs:

* Adapter removal
* Low-quality base removal
* Quality filtering
* Read statistics
* Before/after QC reporting

Trimmed reads are stored under:

```text
results/trimmed/
```

For paired-end data, each sample produces two cleaned FASTQ files:

```text
sample_R1.trimmed.fastq.gz
sample_R2.trimmed.fastq.gz
```

---

# 5. Genome Alignment

Cleaned reads are aligned to the mouse **GRCm39** reference genome using HISAT2.

Conceptually:

```text
Trimmed FASTQ
      │
      ▼
    HISAT2
      │
      ▼
SAM
```

HISAT2 uses the previously generated reference index.

The alignment results are converted into BAM files for downstream processing.

---

# 6. BAM Processing

SAMtools is used to process alignment files.

The workflow performs:

```text
SAM
 ↓
BAM
 ↓
Sorted BAM
 ↓
Indexed BAM
```

Sorted BAM files are stored in:

```text
results/aligned/
```

The BAM index files are generated using:

```text
.bai
```

These files allow efficient access to genomic regions.

---

# 7. Gene Quantification

Gene-level read counting is performed using **featureCounts**.

The annotation file is:

```text
reference/Mus_musculus.GRCm39.116.gtf
```

The pipeline counts reads assigned to genomic features such as:

```text
exon → gene
```

The resulting count files are stored in:

```text
results/counts/
```

The individual sample count files are subsequently combined into a count matrix.

---

# 8. Count Matrix

The count matrix has the following conceptual structure:

```text
Gene       Sample1  Sample2  Sample3  ... Sample30
---------------------------------------------------
Gene_1        120      150      130
Gene_2         45       61       53
Gene_3        900      850      920
...
```

Rows represent genes.

Columns represent biological samples.

The matrix contains **raw integer read counts**.

These raw counts are required by DESeq2.

The matrix should not contain:

* TPM
* FPKM
* CPM
* log-transformed counts
* normalized expression values

---

# 9. Sample Metadata

Metadata are stored in:

```text
metadata/metadata.csv
```

Example:

```text
Sample_Name,cell_line,treatment
CT26_PBS_1,CT26,PBS
CT26_PBS_2,CT26,PBS
CT26_OVM_1,CT26,OVM
CT26_OVM18_1,CT26,OVM18
...
```

The metadata are critical because DESeq2 uses them to determine which samples belong to which experimental groups.

---

# 10. Differential Expression Analysis

Differential expression analysis is performed using **DESeq2**.

The primary design is:

```r
~ cell_line + treatment
```

This design accounts for differences between cell lines while testing the treatment effect.

In simplified form:

```text
Gene expression
       │
       ├── Cell-line effect
       │
       └── Treatment effect
```

This is preferable to ignoring the cell-line factor when multiple cell lines are present.

---

# 11. Treatment Reference

The treatment factor is ordered as:

```r
PBS
OVM
OVM18
```

Therefore:

```text
PBS = reference
```

This allows the primary comparisons to be:

### Comparison 1

```text
OVM vs PBS
```

### Comparison 2

```text
OVM18 vs PBS
```

### Comparison 3

```text
OVM18 vs OVM
```

These comparisons are generated using DESeq2 contrasts.

---

# 12. Low-Count Filtering

Very lowly expressed genes provide limited statistical information.

A typical filtering strategy used in this project is:

```r
keep <- rowSums(counts >= 10) >= 3
```

This keeps genes having at least 10 reads in at least 3 samples.

The purpose is to reduce extremely low-count features before differential expression analysis.

---

# 13. DESeq2 Normalization

DESeq2 performs its own normalization based on sequencing depth and library composition.

The normalized data are used for visualization and interpretation.

Importantly:

```text
Raw counts → DESeq2
```

rather than:

```text
TPM → DESeq2
```

---

# 14. Principal Component Analysis

PCA is performed using variance-stabilized data.

The PCA is used to investigate:

* Sample similarity
* Treatment separation
* Cell-line effects
* Potential outliers
* Batch effects

Example interpretation:

```text
Samples close together
        ↓
Similar global expression profiles
```

Samples that separate strongly may reflect biological differences or technical effects and should be investigated.

---

# 15. Sample Distance Heatmap

A sample-to-sample distance matrix is generated to visualize relationships among the 30 samples.

This can reveal:

* Clustering by treatment
* Clustering by cell line
* Potential outliers
* Unexpected sample relationships

---

# 16. Differential Expression Threshold

For the primary DEG tables, the project uses:

```text
adjusted p-value < 0.05
AND
|log2FoldChange| > 1
```

Therefore:

### Upregulated

```text
padj < 0.05
log2FoldChange > 1
```

### Downregulated

```text
padj < 0.05
log2FoldChange < -1
```

The adjusted p-value is used to account for multiple testing.

---

# 17. Differential Expression Comparisons

The pipeline generates results for:

```text
OVM vs PBS
OVM18 vs PBS
OVM18 vs OVM
```

Each comparison contains:

* Complete DESeq2 results
* Significant genes
* Upregulated genes
* Downregulated genes
* MA plot
* Volcano plot
* Heatmap

---

# 18. MA Plot

MA plots visualize:

```text
log2 Fold Change
        vs
Mean Expression
```

They help identify the relationship between expression level and fold change.

Typical output:

```text
results/deseq2/plots/
```

---

# 19. Volcano Plot

Volcano plots display:

```text
x-axis = log2 Fold Change
y-axis = -log10(adjusted p-value)
```

They provide a visual representation of:

* Strongly upregulated genes
* Strongly downregulated genes
* Statistically significant genes

Highly significant genes with large fold changes appear farther from the center.

---

# 20. Heatmap

Heatmaps are generated for selected differentially expressed genes.

The expression values are typically transformed before visualization.

For example:

```text
Raw counts
    ↓
Variance stabilization
    ↓
Selected DEGs
    ↓
Heatmap
```

Heatmaps help visualize expression patterns across the 30 samples.

---

# 21. Gene Ontology Enrichment

Significant genes are subjected to Gene Ontology enrichment analysis using:

```text
clusterProfiler
```

and:

```text
org.Mm.eg.db
```

GO analysis can investigate biological categories such as:

### Biological Process

Examples include:

```text
immune response
cell proliferation
cell migration
inflammatory response
```

The exact enriched terms depend on the experimental results.

---

# 22. KEGG Pathway Enrichment

KEGG enrichment is performed to identify biological pathways associated with significant genes.

For mouse, the KEGG organism code is:

```text
mmu
```

The analysis is performed separately for:

* Upregulated genes
* Downregulated genes
* Each treatment comparison

---

# 23. Gene ID Conversion

Before GO/KEGG analysis, gene identifiers must be confirmed.

For example, featureCounts may produce:

```text
ENSMUSG00000000001
ENSMUSG00000000028
```

or another identifier type depending on the annotation used.

The correct identifier type must therefore be determined before running:

```r
bitr()
```

For example, if the matrix contains Ensembl IDs:

```r
fromType = "ENSEMBL"
```

If it contains gene symbols:

```r
fromType = "SYMBOL"
```

Incorrect identifier selection can result in poor or empty enrichment results.

---

# 24. Enrichment Background

For rigorous enrichment analysis, the background/universe should ideally represent genes that were actually eligible for DESeq2 testing.

This is preferable to automatically treating the entire mouse genome as the background.

Conceptually:

```text
All annotated mouse genes
        │
        ▼
Genes detected/tested in RNA-seq
        │
        ▼
DESeq2 significant genes
        │
        ▼
GO / KEGG enrichment
```

---

# 25. Expected Results

The final project produces results similar to:

```text
results/
│
├── fastqc/
│
├── multiqc/
│   └── multiqc_report.html
│
├── trimmed/
│
├── aligned/
│   ├── sample.sorted.bam
│   └── sample.sorted.bam.bai
│
├── counts/
│   └── count_matrix.csv
│
├── deseq2/
│   ├── OVM_vs_PBS.csv
│   ├── OVM18_vs_PBS.csv
│   ├── OVM18_vs_OVM.csv
│   │
│   ├── significant/
│   ├── upregulated/
│   ├── downregulated/
│   │
│   └── plots/
│       ├── PCA.png
│       ├── sample_distance.png
│       ├── MA_OVM_vs_PBS.png
│       ├── MA_OVM18_vs_PBS.png
│       ├── MA_OVM18_vs_OVM.png
│       └── volcano_*.png
│
└── enrichment/
    ├── GO/
    └── KEGG/
```

---

# 🚀 Running the Pipeline

## Step 1 — Clone the repository

```bash
git clone <YOUR_GITHUB_REPOSITORY_URL>
cd bulk-rnaseq-snakemake
```

---

## Step 2 — Create the software environment

If an `environment.yml` file is provided:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate bulk-rnaseq
```

---

## Step 3 — Download the data

Run:

```bash
bash scripts/raw_data.sh
```

Verify the FASTQ files:

```bash
ls data/fastq/
```

Verify the reference:

```bash
ls reference/
```

---

# 🧪 Dry Run

Before actually running the pipeline, perform a Snakemake dry run:

```bash
snakemake -n
```

This checks whether Snakemake can determine the required workflow without executing the commands.

For a more detailed dry run:

```bash
snakemake -n -p
```

---

# ▶️ Run the Pipeline

For an 8 GB RAM system, a conservative starting point is:

```bash
snakemake --cores 2
```

To display the commands being executed:

```bash
snakemake --cores 2 -p
```

For systems with more available CPU/RAM:

```bash
snakemake --cores 4
```

The exact number of cores should be adjusted according to available memory.

---

# 🔄 Resume an Interrupted Pipeline

One major advantage of Snakemake is that completed outputs do not normally need to be regenerated.

If the pipeline stops because of an error:

```bash
snakemake --cores 2
```

can be run again after fixing the problem.

Snakemake determines which outputs are already present and which steps still need to be executed.

---

# 🔍 Check the Workflow

To visualize the workflow:

```bash
snakemake --dag | dot -Tpng > workflow.png
```

This produces a directed acyclic graph showing the dependencies between pipeline steps.

---

# 📊 Example Workflow Dependency

```text
FASTQ
 │
 ├── FastQC
 │
 └── fastp
       │
       ▼
   trimmed FASTQ
       │
       ▼
     HISAT2
       │
       ▼
      SAM
       │
       ▼
      BAM
       │
       ▼
  sorted BAM
       │
       ▼
  featureCounts
       │
       ▼
 count matrix
       │
       ▼
    DESeq2
       │
       ├── DEGs
       ├── PCA
       ├── MA plot
       ├── Volcano
       └── Heatmap
             │
             ▼
         GO / KEGG
```

---

# 🧠 Why Snakemake?

Without workflow automation, a 30-sample RNA-seq project can require many manual commands.

For example:

```text
30 samples
×
FastQC
×
fastp
×
HISAT2
×
SAMtools
×
featureCounts
```

Manual execution increases the chance of:

* Typographical errors
* Incorrect sample names
* Missing files
* Inconsistent parameters
* Re-running completed work
* Difficulty reproducing the analysis

Snakemake addresses these problems by defining:

```text
input → processing → output
```

for every analysis step.

---

# ♻️ Reproducibility

The project is designed so that the analysis can be reproduced using:

```text
Snakefile
+
config.yaml
+
metadata
+
scripts
+
software environment
```

The workflow records dependencies between files and automatically determines which steps need to run.

---

# 📌 Important Input Requirements

## Raw reads

The pipeline expects paired-end FASTQ files.

Typical naming:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

---

## Reference genome

Mouse genome:

```text
GRCm39
```

---

## Annotation

The annotation file must correspond appropriately to the reference genome.

Example:

```text
Mus_musculus.GRCm39.116.gtf
```

Genome and annotation versions should be kept consistent.

---

# ⚠️ Large Files Are Not Uploaded to GitHub

The following files should generally **not** be committed to GitHub:

```text
*.fastq
*.fastq.gz
*.sra
*.bam
*.bai
*.sam
*.fa
*.fasta
*.gtf
HISAT2 index files
```

These files can be several GB in size.

Instead, GitHub should contain:

```text
Snakefile
scripts/
metadata/
docs/
README.md
config.yaml
environment.yml
.gitignore
```

---

# 🗂️ `.gitignore`

A suitable `.gitignore` may include:

```gitignore
# Raw sequencing data
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.sra

# Alignment files
*.sam
*.bam
*.bam.bai
*.bai

# Reference genome
*.fa
*.fasta
*.fa.gz
*.fna
*.gtf
*.gtf.gz

# HISAT2 indexes
*.ht2
*.ht2l

# Results
results/

# Temporary files
*.tmp
*.log

# R files
.Rhistory
.RData
.Rproj.user/

# Python
__pycache__/
*.pyc

# OS
.DS_Store
Thumbs.db
```

---

# 📚 Documentation

Additional documentation is provided in the `docs/` directory.

### `REFERENCE.md`

Explains:

* Dataset
* Reference genome
* Annotation
* SRA accessions
* HISAT2 index
* Directory preparation

### `WORKFLOW.md`

Explains:

* Order of pipeline execution
* Snakemake rules
* Input/output relationships
* How to run the workflow

### `DESEQ2_CONNECT.md`

Explains:

* Count matrix
* Metadata
* DESeq2 design
* Contrasts
* DEG filtering
* PCA
* Volcano plots
* Heatmaps
* GO enrichment
* KEGG enrichment

---

# 🧬 Statistical Design

The primary DESeq2 design is:

```r
~ cell_line + treatment
```

This model evaluates treatment differences while accounting for cell-line differences.

The primary treatment comparisons are:

```text
OVM vs PBS
OVM18 vs PBS
OVM18 vs OVM
```

If the scientific question specifically concerns whether treatment effects differ between cell lines, an interaction model can be considered:

```r
~ cell_line * treatment
```

This should only be used when the interaction is part of the biological question.

---

# 📈 Main Questions Addressed

This pipeline can be used to answer questions such as:

### 1. Does OVM change gene expression compared with PBS?

```text
OVM vs PBS
```

### 2. Does OVM18 change gene expression compared with PBS?

```text
OVM18 vs PBS
```

### 3. Does OVM18 differ from OVM?

```text
OVM18 vs OVM
```

### 4. Which genes are significantly upregulated?

```text
padj < 0.05
log2FC > 1
```

### 5. Which genes are significantly downregulated?

```text
padj < 0.05
log2FC < -1
```

### 6. Which biological processes are enriched?

GO enrichment.

### 7. Which biological pathways are affected?

KEGG enrichment.

### 8. Do samples cluster according to treatment or cell line?

PCA and sample-distance analysis.

---

# 🧪 Quality Control Straty

Quality control is performed at multiple stages.

```text
Raw FASTQ
   │
   ▼
FastQC
   │
   ▼
MultiQC
   │
   ▼
fastp
   │
   ▼
FastQC / fastp QC
   │
   ▼
HISAT2 alignment
   │
   ▼
Alignment statistics
   │
   ▼
featureCounts
   │
   ▼
Count matrix QC
   │
   ▼
PCA / sample distance
```

This provides both sequencing-level and expression-level quality assessment.

---

# 🧩 Reproducible Workflow Philosophy

This project follows the principle:

```text
One input
   defined processing step
   ↓
One reproducible output
```

For example:

```text
FASTQ
  ↓
FastQC
  ↓
QC report
```

and:

```text
sorted BAM
  ↓
featureCounts
  ↓
gene counts
```

Snakemake connects these dependencies automatically.

---

# 💻 System Requirements

The pipeline can be run on Linux systems.

Recommended environment:

```text
OS: Ubuntu/Linux
RAM: 8 GB or higher
CPU: Multi-core processor
Storage: Large storage strongly recommended
```

RNA-seq processing generates substantial intees, particularly:

```text
FASTQ
BAM
reference genome
HISAT2 index
```

An 8 GB system can be used for learning and smaller workloads, but processing 30 samples may be slow and storage-intensive.

For production-scale analysis, a workstation, HPC cluster, or cloud environment is preferable.

---

# ⚡ Performance Considerations

For an 8 GB RAM machine, start with:

```bash
snakemake --cores 2
```

Avoid assigning all available CPU threads to memory-intensive processes.

For example, HISAT2 and other alignnt steps can consume substantial memory when multiple samples are processed simultaneously.

Snakemake can be configured to control:

```text
threads
memory
parallel jobs
```

---

# 🧪 Troubleshooting

## Check Snakemake installation

```bash
snakemake --version
```

## Check HISAT2

```bash
hisat2 --version
```

## Check featureCounts

```bash
featureCounts -v
```

## Check SAMtools

```bash
samtools --version
```

## Check R

```bash
R --version
```

## Check DESeq2

Inside R:

```r
library(DESeq2)
```## Check clusterProfiler

```r
library(clusterProfiler)
```

## Check mouse annotation database

```r
library(org.Mm.eg.db)
```

---

# 🔬 Important Reproducibility Checks

Before running DESeq2:

```r
all(rownames(metadata) == colnames(counts))
```

This should return:

```text
TRUE
```

Also verify:

```r
all(counts == round(counts))
```

This should return:

```text
TRUE
```

The count matrix should contain raw integer counts.

---

# 📊 Final Deliverables

The project generates the following major ots:

### QC

```text
FastQC reports
MultiQC report
```

### Alignment

```text
Sorted BAM
BAM index
```

### Quantification

```text
Gene-level count matrix
```

### Differential expression

```text
DESeq2 result tables
Significant DEG tables
Upregulated genes
Downregulated genes
```

### Visualization

```text
PCA
Sample-distance heatmap
MA plots
Volcano plots
DEG heatmaps
```

### Functional analysis

```text
GO enrichment
KEGG enrichment
```

---

# 🎯 Portfolio Objectives

This project demonstrates prtical experience with:

* Bulk RNA-seq analysis
* Linux
* Bash scripting
* Snakemake
* Workflow automation
* SRA data handling
* FASTQ processing
* Quality control
* Read trimming
* HISAT2 alignment
* SAMtools
* featureCounts
* R
* DESeq2
* Statistical analysis
* Differential expression
* Gene Ontology
* KEGG pathway analysis
* Data visualization
* Reproducible research

---

# 💼 Skills Demonstrated

This project demonstrates the ability to work with a workflow such as:

```text
Public sequencing data
      ↓
Linux/Bash
        ↓
Workflow automation
        ↓
NGS quality control
        ↓
RNA-seq alignment
        ↓
Read quantification
        ↓
Statistical analysis
        ↓
Biological interpretation
```

It therefore represents an end-to-end bulk RNA-seq analysis rather than an isolated command-line exercise.

---

# 📌 Project Limitations

This project has several important limitations.

1. Differential expression results depend on the quality and biological design of the original dataset.

2. The selected thresholds:

```text
padj < 0.05
|log2FC| > 1
```

are analysis choices and should not be interpreted as universal biological cutoffs.

3. GO/KEGG results depend on correct gene identifier conversion and appropriate background selection.

4. The additive model:

```r
~ cell_line + treatment
```

does not directly test treatment-by-cell-line interaction.

5. Alignment and quantification results depend on the selected reference genome and annotation versions.

6. A successful computational pipeline does not by itself establish biological causality; experimental validation may be required.

---

# 🔁 Reproducing the Project

The basic workflow is:

```bash
git clone <YOUR_GITHUB_REPOSITORY_URL>

cd bulk-rnaseq-snakemake

conda env create -f environment.yml

conda activate bulk-rnaseq

bash scripts/raw_data.sh

snakemake -n

snakemake --cores 2
```

The exact commands may vary depending on the environment and configuration.

---

# 👨‍💻 Author

**Parth Vyas**

MSc Bioinformatics

Thisct was developed as a practical demonstration of an automated, reproducible bulk RNA-seq workflow using Snakemake and R.

---

# 📜 License

This project is intended for educational and research purposes.

Add an appropriate open-source license to the repository if you plan to distribute the code publicly.

---

# ⭐ Acknowledgements

This project makes use of publicly available sequencing data and open-source bioinformatics software.

Major tools used include:

* Snakemake
* SRA Toolkit
* FastQC
* Multi fastp
* HISAT2
* SAMtools
* featureCounts
* DESeq2
* clusterProfiler
* R

---

# 📌 Quick Summary

```text
Dataset
    GSE334005

Organism
    Mus musculus

Genome
    GRCm39

Samples
    30

Design
    3 treatments × 3 cell lines

Treatments
    PBS / OVM / OVM18

Read type
    Paired-end

Workflow
    SRA
      ↓
    FASTQ
      ↓
    FastQC
      ↓
    MultiQC
      ↓
    fastp
      ↓
    HISAT2
      ↓
    SAMtools
      ↓
    featureCounts
      ↓
    Count Matrix
      ↓
    DES     ↓
    DEG Analysis
      ↓
    GO
      ↓
    KEGG
```

**This repository provides a reproducible end-to-end framework for processing 30 mouse bulk RNA-seq samples from raw sequencing data through differential expression and functional enrichment analysis.**
