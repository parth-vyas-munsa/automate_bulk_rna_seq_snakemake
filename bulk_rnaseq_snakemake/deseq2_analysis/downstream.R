library(DESeq2)

counts = read.delim("GSE334005_RNAseq_CT26_B16F10_LLC_read.tsv")
##write.csv(counts , "counts.csv" , row.names = FALSE )
head(counts)

metadata = read.csv("SraRunTable (2).csv")
##write.csv(metadata , "matadata.csv" , row.names = FALSE )
head(metadata)
genes = counts[, c("gene_id" , "gene_symbol")]
counts = counts[,-c(1, 2)]
head(counts)

rownames(metadata) <- metadata$Sample_Name
head(metadata)

metadata = metadata[,c("treatment") , drop = FALSE]
table(metadata$treatment)
metadata


dim(metadata)
colnames(counts)
colnames(counts)
rownames(metadata)
all(colnames(counts) %in% rownames(metadata))

metadata <- metadata[colnames(counts), ]
all(rownames(metadata) == colnames(counts))
head(metadata)

table(metadata$cell_line)
table(metadata$cell_line, metadata$treatment)

metadata$cell_line <- factor(
  metadata$cell_line
)

metadata$treatment <- factor(
  metadata$treatment,
  levels = c(
    "PBS",
    "OVM (M1)",
    "OVM18 (M1-mutIL18)"
  )
)

levels(metadata$treatment)

is.numeric(as.matrix(counts))
any(counts < 0)
all(counts == round(counts))
##filtering low counts
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]
dim(counts)
dim(counts_filtered)
any(counts_filtered %% 1 != 0)
counts_filtered <- round(counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ cell_line + treatment
)
dds
design(dds)

dds <- DESeq(dds)
ds_resultsNames(dds)

res_OVM_vs_PBS <- ds_results(
  dds,
  contrast = c(
    "treatment",
    "OVM (M1)",
    "PBS"
  ),
  alpha = 0.05
)
head(res_OVM_vs_PBS)

res_OVM18_vs_PBS <- ds_results(
  dds,
  contrast = c(
    "treatment",
    "OVM18 (M1-mutIL18)",
    "PBS"
  ),
  alpha = 0.05
)
head(res_OVM18_vs_PBS)
write.csv(res_OVM18_vs_PBS , "res_OVM18_vs_PBS_dds" , row.names = FALSE)


res_OVM18_vs_OVM <- ds_results(
  dds,
  contrast = c(
    "treatment",
    "OVM18 (M1-mutIL18)",
    "OVM (M1)"
  ),
  alpha = 0.05
)
head(res_OVM18_vs_OVM)
write.csv(res_OVM18_vs_OVM , "res_OVM18_vs_OVM_dds.csv" , row.names = FALSE)

##ovm vs pbs
res_OVM_vs_PBS <- res_OVM_vs_PBS[
  order(res_OVM_vs_PBS$padj),
]
head(res_OVM_vs_PBS)


res_OVM_vs_PBS_clean <- res_OVM_vs_PBS[
    !is.na(res_OVM_vs_PBS$padj),
]

##ovm vs ovm 18
res_OVM18_vs_OVM <- res_OVM18_vs_OVM[
  order(res_OVM18_vs_OVM$padj),
]
head(res_OVM18_vs_OVM)


res_OVM18_vs_OVM_clean <- res_OVM18_vs_OVM[
  !is.na(res_OVM18_vs_OVM$padj),
]

##ovm 18 vs pbs
res_OVM18_vs_PBS <- res_OVM18_vs_PBS[
  order(res_OVM18_vs_PBS$padj),
]
head(res_OVM18_vs_PBS)


res_OVM18_vs_PBS_clean <- res_OVM18_vs_PBS[
  !is.na(res_OVM18_vs_OVM$padj),
]
head(res_OVM18_vs_PBS_clean)

##Define significant DEGs
sig_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
nrow(sig_OVM_vs_PBS)
##Upregulated
up_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM_vs_PBS)
nrow(down_OVM_vs_PBS)


##OVM18_vs_PBS
sig_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
##Upregulated
up_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM18_vs_PBS)
nrow(down_OVM18_vs_PBS)

##_OVM18_vs_OVM

sig_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

##Upregulated
up_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM18_vs_OVM)
nrow(down_OVM18_vs_OVM)


dir.create("ds_results", showWarnings = FALSE)

write.csv(
  as.data.frame(res_OVM_vs_PBS),
  "ds_results/DESeq2_OVM_vs_PBS.csv"
)

write.csv(
  as.data.frame(res_OVM18_vs_PBS),
  "ds_results/DESeq2_OVM18_vs_PBS.csv"
)

write.csv(
  as.data.frame(res_OVM18_vs_OVM),
  "ds_results/DESeq2_OVM18_vs_OVM.csv"
)

write.csv(
  as.data.frame(sig_OVM_vs_PBS),
  "ds_results/significant_OVM_vs_PBS.csv"
)

write.csv(
  as.data.frame(sig_OVM18_vs_PBS),
  "ds_results/significant_OVM18_vs_PBS.csv"
)

write.csv(
  as.data.frame(sig_OVM18_vs_OVM),
  "ds_results/significant_OVM18_vs_OVM.csv"
)



#######PCA



vsd <- vst(
  dds,
  blind = FALSE
)


plotPCA(
  vsd,
  intgroup = c("treatment", "cell_line")
)

png(
  "ds_results/PCA.png",
  width = 1000,
  height = 800
)

plotPCA(
  vsd,
  intgroup = c("treatment", "cell_line")
)

dev.off()
BiocManager::install('pheatmap')
library(pheatmap)

sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)
pheatmap(
  sample_dist_matrix,
  annotation_col = metadata[, c(
    "cell_line",
    "treatment"
  )],
  filename = "ds_results/sample_distance_heatmap.png"
)

plot(
  hclust(sample_dist)
)
png(
  "ds_results/sample_clustering.png",
  width = 1000,
  height = 800
)

plot(
  hclust(sample_dist),
  main = "Sample clustering"
)

dev.off()
png(
  "ds_results/MA_OVM_vs_PBS.png",
  width = 1000,
  height = 800
)

plotMA(
  res_OVM_vs_PBS,
  ylim = c(-5, 5)
)

dev.off()

#########volcano plot for res_OVM_vs_PBS
volcano_data <- as.data.frame(res_OVM_vs_PBS)

volcano_data$significance <- "Not significant"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 1
] <- "Up"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < -1
] <- "Down"

library(ggplot2)

ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "OVM vs PBS",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )

volcano_data$gene <- rownames(volcano_data)
top_genes <- head(
  volcano_data[
    order(volcano_data$padj),
  ],
  10
)
library(ggrepel)
ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene)
  ) +
  theme_classic()

top50 <- head(
  rownames(
    res_OVM_vs_PBS[
      order(res_OVM_vs_PBS$padj),
    ]
  ),
  50
)

mat <- assay(vsd)[top50, ]
mat_scaled <- t(
  scale(t(mat))
)
annotation <- metadata[, c(
  "cell_line",
  "treatment"
)]
pheatmap(
  mat_scaled,
  annotation_col = annotation,
  show_rownames = TRUE,
  fontsize_row = 6,
  filename = "ds_results/heatmap_top50_OVM_vs_PBS.png"
)
normalized_counts <- counts(
  dds,
  normalized = TRUE
)
write.csv(
  as.data.frame(normalized_counts),
  "ds_results/normalized_counts.csv"
)


#########volcano plot for res_OVM_vs_OVM18
library(ggrepel)
library(ggplot2)
volcano_data <- as.data.frame(res_OVM18_vs_OVM)

volcano_data$significance <- "Not significant"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 1
] <- "Up"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < -1
] <- "Down"

ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "OVM vs OVM18",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )
volcano_data$gene <- rownames(volcano_data)


top_genes <- head(
  volcano_data[
    order(volcano_data$padj),
  ],
  10
)

ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene)
  ) +
  theme_classic()
top50 <- head(
  rownames(
    res_OVM18_vs_OVM[
      order(res_OVM18_vs_OVM$padj),
    ]
  ),
  50
)

mat <- assay(vsd)[top50, ]

annotation <- metadata[, c(
  "cell_line",
  "treatment"
)]

pheatmap(
  mat_scaled,
  annotation_col = annotation,
  show_rownames = TRUE,
  fontsize_row = 6,
  filename = "ds_results/heatmap_top50_OVM_vs_OVM18.png"
)

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  as.data.frame(normalized_counts),
  "ds_results/normalized_counts_ovm_vs_ovm18.csv"
)

png(
  "ds_results/dispersion_plot_ovm_vs_ovm18.png",
  width = 1000,
  height = 800
)

plotDispEsts(dds)

dev.off()

#########volcano plot for res_PBS_vs_OVM18
library(ggrepel)
library(ggplot2)
volcano_data <- as.data.frame(res_OVM18_vs_PBS)

volcano_data$significance <- "Not significant"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 1
] <- "Up"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < -1
] <- "Down"


ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    title = "OVM_18 vs PBS",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  )

volcano_data$gene <- rownames(volcano_data)

top_genes <- head(
  volcano_data[
    order(volcano_data$padj),
  ],
  10
)
ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  geom_point(alpha = 0.6) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene)
  ) +
  theme_classic()

top50 <- head(
  rownames(
    res_OVM_vs_PBS[
      order(res_OVM18_vs_PBS$padj),
    ]
  ),
  50
)
mat <- assay(vsd)[top50, ]

mat_scaled <- t(
  scale(t(mat))
)

annotation <- metadata[, c(
  "cell_line",
  "treatment"
)]

pheatmap(
  mat_scaled,
  annotation_col = annotation,
  show_rownames = TRUE,
  fontsize_row = 6,
  filename = "ds_results/heatmap_top50_OVM18_vs_PBS.png"
)

normalized_counts <- counts(
  dds,
  normalized = TRUE
)


write.csv(
  as.data.frame(normalized_counts),
  "ds_results/normalized_counts_OVM_18_VS_PBS.csv"
)





BiocManager::install("org.Mm.eg.db")

############go enrichment
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)


##### go enrichment of all samples  
#####go enrich ment for ovm_vs_pbs ##########################################

sig_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
nrow(sig_OVM_vs_PBS)
##Upregulated
up_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM_vs_PBS)
nrow(down_OVM_vs_PBS)


up_genes <- rownames(up_OVM_vs_PBS)
down_genes <- rownames(down_OVM_vs_PBS)

head(up_genes)
head(down_genes)

ego_up <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_up)
write.csv(ego_up , "ds_results/ego_uo_OVM_VS_PBS.csv")

ego_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_down)
write.csv(ego_up , "ds_results/ego_down_OVM_VS_PBS.csv")


##### go enrich ment for ovm18 vs pbs


##OVM18_vs_PBS
sig_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
##Upregulated
up_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM18_vs_PBS)
nrow(down_OVM18_vs_PBS)

up_genes <- rownames(up_OVM18_vs_PBS)
down_genes <- rownames(down_OVM18_vs_PBS)

head(up_genes)
head(down_genes)



ego_up <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_up)
write.csv(ego_up , "ds_results/ego_uo_OVM18_VS_PBS.csv")

ego_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_up)
write.csv(ego_up , "ds_results/ego_down_OVM18_VS_PBS.csv")









##_OVM18_vs_OVM go enrich ment

sig_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

##Upregulated
up_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM18_vs_OVM)
nrow(down_OVM18_vs_OVM)

up_genes <- rownames(up_OVM18_vs_OVM)
down_genes <- rownames(down_OVM18_vs_OVM)

head(up_genes)
head(down_genes)



ego_up <- enrichGO(
  gene          = up_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_up)
write.csv(ego_up , "ds_results/ego_uP_OVM18_VS_OVM.csv")

ego_down <- enrichGO(
  gene          = down_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_down)
write.csv(ego_down , "ds_results/ego_down_OVM18_VS_OVM.csv")


dds_LRT <- DESeq(
  dds,
  test = "LRT",
  reduced = ~ cell_line
)

res_LRT <- ds_results(dds_LRT)

res_LRT <- res_LRT[
  order(res_LRT$padj),
]





###kegg pathway for all samples

##OVM18_vs_PBS
sig_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
##Upregulated
up_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM18_vs_PBS <- subset(
  res_OVM18_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM18_vs_PBS)
nrow(down_OVM18_vs_PBS)

up_genes <- rownames(up_OVM18_vs_PBS)
down_genes <- rownames(down_OVM18_vs_PBS)

kegg_up <- enrichKEGG(
  gene = up_genes,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(kegg_up)
write.csv(kegg_up , "results/kegup_ovm18_vs_pbs.csv")

kegg_down <- enrichKEGG(
  gene = down_genes,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(kegg_down)
write.csv(kegg_down , "results/keggdown_ovm18_vs_pbs.csv")




###########kegg for ovm vs ovm 18

sig_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

##Upregulated
up_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM18_vs_OVM <- subset(
  res_OVM18_vs_OVM,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM18_vs_OVM)
nrow(down_OVM18_vs_OVM)

up_genes <- rownames(up_OVM18_vs_OVM)
down_genes <- rownames(down_OVM18_vs_OVM)


kegg_up <- enrichKEGG(
  gene = up_genes,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(kegg_up)
write.csv(kegg_up , "results/kegup_ovm18_vs_ovm.csv")

kegg_down <- enrichKEGG(
  gene = down_genes,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(kegg_down)
write.csv(kegg_down , "results/keggdown_ovm18_vs_ovm.csv")







######kegg for ovm vs pbs

sig_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
nrow(sig_OVM_vs_PBS)
##Upregulated
up_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)
##Downregulated
down_OVM_vs_PBS <- subset(
  res_OVM_vs_PBS,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange < -1
)
nrow(up_OVM_vs_PBS)
nrow(down_OVM_vs_PBS)


up_genes <- rownames(up_OVM_vs_PBS)
down_genes <- rownames(down_OVM_vs_PBS)


kegg_up <- enrichKEGG(
  gene = up_genes,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(kegg_up)
write.csv(kegg_up , "results/kegup_ovm_vs_pbs.csv")












