# DESeq2, GO and KEGG Analysis

## Project

This document describes the downstream analysis of the 30-sample mouse bulk RNA-seq dataset.

The workflow starts after `featureCounts` and ends with:

- Count matrix
- DESeq2 differential expression
- Significant genes
- PCA
- Sample clustering
- MA plots
- Volcano plots
- Heatmaps
- GO enrichment
- KEGG enrichment
- Final result tables

---

# 1. Analysis Overview

The downstream workflow is:

```text
featureCounts
     │
     ▼
Individual count files
     │
     ▼
Combine count files
     │
     ▼
Count matrix
     │
     ├──────────────┐
     │              │
     ▼              ▼
Metadata        Count QC
     │              │
     └──────┬───────┘
            ▼
          DESeq2
            │
            ▼
     Differential genes
            │
      ┌─────┴─────┐
      ▼           ▼
     Up          Down
      │           │
      └─────┬─────┘
            ▼
       GO enrichment
            │
            ▼
      KEGG enrichment


the deg perform between 3groups equally separate 
OVM VS OVM18 CELL LINE
OVM VS PBS CELL LINE
PBS VS OVM18 CELLLINE

