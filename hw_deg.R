# ==============================================================================
# Many DEG related packages come from BioConductor
# ==============================================================================
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


# Install DESeq2 and GO packages.
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")

install.packages("pheatmap")
install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("cowplot")

# ==============================================================================
# Load required libraries.
# ==============================================================================
library(readr)
library(tidyr)
library(dplyr)


# ==============================================================================
# Read TSV dataset from the file
# And convert it into a matrix (by removing the first column "ENSEMBLID").
# Then replace row indexes with Ensembl IDs.
# ==============================================================================
counts = read_tsv('~/hw/counts/counts_data.tsv')
counts_mat = counts %>% 
  select(-ENSEMBLID) %>% 
  as.matrix()
# trim the version suffix from the gene id, otherwise GSEA will not work.
rownames(counts_mat) <- lapply(counts$ENSEMBLID, sub, pattern = "\\.\\d+$", replacement = "")
# rownames(counts_mat) <- counts$ENSEMBLID
mode(counts_mat) <- "integer"

# ==============================================================================
# Prepare metadata and it's data-frame.
# ==============================================================================
meta_data = tribble(
  ~sample, ~senescence,
  'SRR1660555', 'Young',
  'SRR1660556', 'Young',
  'SRR1660557', 'Young',
  'SRR1660558', 'Old',
  'SRR1660559', 'Old',
  'SRR1660560', 'Old'
)
meta_data_df = meta_data %>% 
  select(-sample, ) %>%
  as.data.frame()
rownames(meta_data_df) = meta_data$sample
meta_data_df$senescence <- factor(meta_data_df$senescence, levels = c('Young', 'Old'))

# ==============================================================================
# Compute the pairwise Pearson correlation between all variables (reads).
# E.g. we get the covariance matrix (between samples).
#
# Simultaneously we scale the data, to account for differences in
# gene expression levels.
# ==============================================================================
cors = cor(log(1 + counts_mat), m='p')

# ==============================================================================
# Draw the heat-map
# ==============================================================================
library(pheatmap)
pheatmap(cors, 
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         annotation_col = meta_data_df,
         annotation_colors = list(senescence = c('Young' = "green", 'Old' = "yellow")),
         main = 'Pearson')

# On the heat-map we can clearly see that Old and Young samples group together.
# Samples within a group a red-ish, what means they correlate
# and samples between groups are anticorrelate.

# TODO: uncomment to export the image
# png("cors.png", width = 800, height = 600)
# dev.off()


# ==============================================================================
# DEseq2
# ==============================================================================
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = meta_data_df,
  design = ~senescence)

# gene-dependent normalization factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds) %>% as.numeric()

# And plot PCA.
# Here I am using vst to account for a variance that I may observe because
# a gene is more expressed and therefore has a greater variance.
plotPCA(vst(dds), intgroup = "senescence")

# As you can see groups were clearly formed.

# Perform acutall analysis.
dds <- DESeq(dds)
# lists the coeficients
resultsNames(dds)
# and extract the result table
res = results(dds, name = "senescence_Old_vs_Young", alpha = 0.1)
summary(res)

# Now we can plot the Volcano chart, that visualizes which genes conform to
# log2FoldChange and p-value criteria (in red).
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

# 1. convert the object into data frame
# 2. then add the gene column
# 3. keep only gene, log2FoldChnage, padj columns
# 4. drop all padj N/As
# 5. and keep only if adjusted p-value is good enough for us
# 6. sort by the absolute value of log2FC
res2 <- res %>% as.data.frame() %>%
  mutate(gene = rownames(.)) %>% 
  select(gene, log2FoldChange, padj) %>% 
  drop_na(padj) %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(abs(log2FoldChange))) 


# Or I can plot counts for the best and the worst gene from the dataset
# filtered by the adjusted p-value.
plotCounts(dds, gene=head(res2, n=1)[[1]], intgroup="senescence")
plotCounts(dds, gene=tail(res2, n=1)[[1]], intgroup="senescence")

# The log2FC is clearly smaller for the latest gene.


# ==============================================================================
# Over Representation Analysis
#
# Get UP and DOWN regulated genes separately.
# ==============================================================================
library(org.Hs.eg.db)

res_up <- res2 %>%
  filter(log2FoldChange > 1)  %>%
  top_n(1000)
nrow(res_up)

res_down <- res2 %>%
  filter(log2FoldChange < -1)  %>%
  top_n(1000)
nrow(res_down)

ego_up <- clusterProfiler::enrichGO(
  gene = res_up$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  minGSSize = 50,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  # Use full TERM2GENE table
  # universe = res$gene,
  readable = TRUE)
nrow(ego_up)

ego_down <- clusterProfiler::enrichGO(
  gene = res_down$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  minGSSize = 50,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE)
nrow(ego_down)

library(cowplot)
p1 <- clusterProfiler::dotplot(ego_up, showCategory = 10) + ggplot2::ggtitle("Up Regulated")
p2 <- clusterProfiler::dotplot(ego_down, showCategory = 10) + ggplot2::ggtitle("Down Regulated")
cowplot::plot_grid(p1, p2)


clusterProfiler::goplot(ego_up, showCategory = 10)


edox <- clusterProfiler::setReadable(ego, "org.Hs.eg.db", "ENSEMBL")
p1 <- clusterProfiler::cnetplot(edox, categorySize="pvalue", cex_label_gene = 0.8)
p2 <- clusterProfiler::heatplot(edox, showCategory=5)
cowplot::plot_grid(p1, p2)


# ==============================================================================
# Gene Set Enrichment Analysis 
#
# ORA may fail where the difference in
# gene expression is small, but evidenced in coordinated way in a set of
# related genes.
#
# GSEA aggregates the per gene statistics across genes within a gene set,
# therefore making it possible to detect situations where all genes in a
# predefined set change in a small but coordinated way.
# ==============================================================================


# GSEA requires specific format of the input geneList
# See [https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist] for more details.

res3 <- res_up

geneIds <- clusterProfiler::bitr(
  res3$gene,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  drop = FALSE,
  OrgDb = org.Hs.eg.db)$ENTREZID
geneIds[duplicated(geneIds)] <- NA
geneList <- res3$log2FoldChange[!is.na(geneIds)]
names(geneList) <- geneIds[!is.na(geneIds)]
geneList <- sort(geneList, decreasing = TRUE)

ego_up <- clusterProfiler::gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP", # CC MF
  # exponent = 0.2,
  keyType = "ENTREZID", # "ENSEMBL",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  # pAdjustMethod = "BH",
  # seed = TRUE,
  # by = "fgsea",
  verbose = TRUE)
nrow(ego_up)


res3 <- res %>% as.data.frame() %>%
  mutate(gene = rownames(.)) %>% 
  dplyr::select(gene, log2FoldChange, padj, lfcSE) %>% 
  drop_na(padj, lfcSE) %>% 
  filter(padj < 0.05) %>% 
  filter(log2FoldChange < -2) %>%
  arrange(desc(abs(log2FoldChange)))

geneIds <- clusterProfiler::bitr(
  res3$gene,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  drop = FALSE,
  OrgDb = org.Hs.eg.db)$ENTREZID
geneIds[duplicated(geneIds)] <- NA
geneList <- abs(res3$log2FoldChange[!is.na(geneIds)])
names(geneList) <- geneIds[!is.na(geneIds)]
geneList <- sort(geneList, decreasing = TRUE)
length(geneList)
ego_down <- clusterProfiler::gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "BP", # CC MF
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE)
nrow(ego_down)


clusterProfiler::dotplot(ego_up, showCategory=10) + ggplot2::ggtitle("Up regulation")
clusterProfiler::dotplot(ego_down, showCategory=10) + ggplot2::ggtitle("Down regulation")


# ==============================================================================
# KEGG Enrichment
# ==============================================================================

gene <- res2$gene[abs(res2$log2FoldChange) > 2]
gene <- clusterProfiler::bitr(
  gene,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  drop = FALSE,
  OrgDb = org.Hs.eg.db)$ENTREZID

kk <- clusterProfiler::enrichKEGG(
  gene = gene,
  organism     = 'hsa',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500)
nrow(kk)

clusterProfiler::dotplot(kk, showCategory = 10)

# ==============================================================================
# Export the list of genes to the text file.
# ==============================================================================
gene_ids <- clusterProfiler::bitr(
  res2$gene,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  drop = FALSE,
  OrgDb = org.Hs.eg.db)$SYMBOL

lapply(gene_ids[!is.na(gene_ids)], write, "gene-list.txt", append=TRUE, ncolumns=1000)

# library(pheatmap)
# library(DOSE)
# library(pathview)
# library(ggplot2)