library(dplyr)
library(tidyverse)
library(tibble)
library(tidyr)
library(Seurat)
library(viridis)
library(RColorBrewer)
library(sceasy)
library(escape)
library(openxlsx)

# Load 'Esteso et al' BCG Seurat R file and subset Oncotice-treated cells (GSE203098)
bcg <- readRDS("D:/sc/blood/bcg_esteso.rds")
Idents(bcg) <- 'condition'
oncotice <- subset(bcg, idents=c('Oncotice'))
n <- 36
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# FIGURE 2A
DimPlot(oncotice, label=TRUE, cols=col_vector)
# Subset oncotice-treated nk cells
Idents(oncotice) <- 'annotation'
new.cluster.ids <- c("Rest", "Rest", "Rest", "NK cell", "Rest", "Rest", "Rest", "Rest", "Rest", "Rest", "Rest", "Rest", "Rest", "Rest")
names(new.cluster.ids) <- levels(oncotice)
oncotice <- RenameIdents(oncotice, new.cluster.ids)
# FIGURE 2B
VlnPlot(oncotice, features=c('TYROBP', 'NCAM1', 'CD3D'), cols=c('gray40','#FFFF99'))
# Clustering and annotation
nk <- subset(oncotice, idents=c('NK cell'))
nk <- NormalizeData(nk)
nk <- FindVariableFeatures(nk, selection.method = "mean.var.plot")
nk <- ScaleData(nk, features = rownames(nk), block.size = 100) 		
nk <- RunPCA(nk, features = VariableFeatures(object = nk))
ElbowPlot(nk)
nk <- FindNeighbors(nk, dims = 1:15)
nk <- FindClusters(nk, resolution = 0.5)
nk <- RunUMAP(nk, dims = 1:15, umap.method = 'umap-learn', metric = 'correlation', spread=1, min.dist=0.05)
new.cluster.ids <- c("1", "2", "3")
names(new.cluster.ids) <- levels(nk)
nk <- RenameIdents(nk, new.cluster.ids)
number <- Idents(nk)
number <- unname(number)
nk@meta.data$number = number

# Load unstimulated PBMCs at day 0 from Lee JS et al (GSE149689)
pbmc <- Read10X(data.dir = "D:/sc/blood/control/raw/filtered/")
pbmc <- CreateSeuratObject(counts = pbmc)
metadata <- read.table('meta.tsv', sep='\t', header=TRUE)	## load metadata file from original paper
cells.use <- metadata$cellId
pbmc <- subset(pbmc, cells = cells.use)
pbmc@meta.data <- metadata
Idents(pbmc) <- 'Disease.group'
pbmc <- subset(pbmc, idents=c('Healthy Donor'))
embeds <- read.table('tSNE.coords.tsv.gz', sep='\t')	## load tSNE coordinates from original paper
row.names(embeds) <- embeds$V1
embeds <- embeds[,-1]
common <- intersect(rownames(embeds), colnames(pbmc))
embeds <- embeds[match(common, rownames(embeds)),]
identical(colnames(pbmc), rownames(embeds))
metadata <- read.table('meta.tsv', sep='\t', header=TRUE)
metadata2 <- pbmc@meta.data
metadata2$cellId <- colnames(pbmc)
common <- intersect(metadata$cellId, metadata2$cellId)
metadata <- metadata[match(common, metadata$cellId),]
identical(metadata$cellId, metadata2$cellId)
pbmc@meta.data <- metadata
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "mean.var.plot")
pbmc <- ScaleData(pbmc, features = rownames(pbmc), block.size = 100) 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunTSNE(pbmc, dims = 1:30)
embeds <- embeds %>% rename(tSNE_1 = V2, tSNE_2 = V3)	
embeds <- data.matrix(embeds)
pbmc[["tsne"]]@cell.embeddings <- embeds
highlight <- WhichCells(pbmc, idents=c("NK cell"))
# Supplementary Figure 2B
DimPlot(pbmc, cells.highlight= highlight, cols.highlight = c("grey60"), cols= "lightgrey", sizes.highlight = 0.3)
# Subset unstimulated NK cells at day 0
nk_control <- subset(pbmc, idents=c('NK cell'))
metadata <- read.table('meta.tsv', sep='\t', header=TRUE)
metadata2 <- nk_control@meta.data
metadata2$cellId <- colnames(nk_control)
common <- intersect(metadata$cellId, metadata2$cellId)
metadata <- metadata[match(common, metadata$cellId),]
identical(metadata$cellId, metadata2$cellId)
metadata <- metadata %>% rename(percent.mt = Percentage.of.mitochondrial.gene, nFeature_RNA = Number.of.Gene, nCount_RNA = Number.of.UMI)
nk_control@meta.data <- metadata

# Merge Oncotice-stimulated nk cells and unstimulated nk cells
metadata <- nk@meta.data
metadata <- metadata %>% mutate(Sample.ID = case_when(endsWith(HTO_classification, "C0251-TotalSeqC") ~ "BCG 1", endsWith(HTO_classification, "C0252-TotalSeqC") ~ "BCG 2", endsWith(HTO_classification, "C0253-TotalSeqC") ~"BCG 3"))
nk@meta.data <- metadata
nk.merged <- merge(nk, y = c(nk_control), add.cell.ids = c("bcg", "control"), project = "blood")
features_keep <- row.names(nk_control)
features_keep <- as.vector(features_keep)
nk.merged <- subset(nk.merged, features = features_keep)
# Clustering and annotation
nk.merged <- NormalizeData(nk.merged)
nk.merged <- FindVariableFeatures(nk.merged, selection.method = "mean.var.plot")
nk.merged <- ScaleData(nk.merged, vars.to.regress=c('Sample.ID','percent.mt'), features = rownames(nk.merged))
nk.merged <- RunPCA(nk.merged, features = VariableFeatures(object = nk.merged))
nk.merged <- FindNeighbors(nk.merged, dims = 1:15)
nk.merged <- FindClusters(nk.merged, resolution = 0.2)
nk.merged <- RunUMAP(nk.merged, dims = 1:15, spread=1, min.dist=0.05)
DimPlot(nk.merged, reduction = "umap", label=TRUE)
new.cluster.ids <- c("4", "5", "1", "6", "2")
names(new.cluster.ids) <- levels(nk.merged)
nk.merged <- RenameIdents(nk.merged, new.cluster.ids)
cluster <- Idents(nk.merged)
cluster <- unname(cluster)
nk.merged@meta.data$cluster = cluster
metadata <- nk.merged@meta.data
metadata$cluster <- coalesce(metadata$number, metadata$cluster)
nk.merged@meta.data <- metadata
Idents(nk.merged) <- 'cluster'
my_levels <- c("1", "2", "3", "4", "5", "6")
Idents(nk.merged) <- factor(Idents(nk.merged), levels= my_levels)
cluster <- Idents(nk.merged)
cluster <- unname(cluster)
nk.merged@meta.data$cluster = cluster
# FIGURE 2C
DimPlot(nk.merged, cols=c('#457b9d', '#ffb942', '#e63946', 'gray40', 'gray60', 'gray80'))

# FIGURE 2D
DoHeatmap(nk, features = c('CCL5', 'NKG7', 'TYROBP', 'KLRB1', 'CST7', 'MKI67', 'TOP2A', 'RASGRP2', 'KLF2', 'IGFBP7', 'HSP90AB1', 'FABP5', 'BATF3', 'ZBED2', 'TNFRSF4'), size = 2, group.bar.height = 0.02, group.colors = c('#457b9d', '#ffb942', '#e63946')) + NoLegend() + scale_fill_gradientn(colors = c('white',"white", "black"))

# FIGURE 2E
cluster1.markers <- FindMarkers(nk, ident.1 = 1, min.pct = 0.25)
cluster1 <- filter(cluster1.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster2.markers <- FindMarkers(nk, ident.1 = 2, min.pct = 0.25)
cluster2 <- filter(cluster2.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster3.markers <- FindMarkers(nk, ident.1 = 3, min.pct = 0.25)
cluster3 <- filter(cluster3.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
metascape_input <- data.frame(C1=rownames(cluster1), C2=rownames(cluster2), C3=rownames(cluster3))
write.xlsx(metascape_input, 'metascape_input.xlsx')	## run metascape (https://metascape.org/gp/index.html#/main/step1)
GO_dotplot <- read.csv('GO_dotplot.csv', header=TRUE, sep=';')	## load in-house csv generated from metascape
library(ggplot2)
GO_dotplot$Cluster <- as.character(GO_dotplot$Cluster)
GO_dotplot$GeneRatio[GO_dotplot$GeneRatio == 0] <- NA
GO_dotplot$logp.adjust[GO_dotplot$logp.adjust == 0] <- NA
ggplot(data = GO_dotplot, aes(x = Cluster, y = GO)) + geom_point(aes(fill=-(GeneRatio), size=-(logp.adjust)), colour="black", shape=21, stroke=0.5) + theme_bw() + ylab("") + xlab("Clusters") + ggtitle("GO enrichment analysis") + scale_fill_continuous(high = "#132B43", low = "#56B1F7", limits = c(-20, 0)) 

# FIGURE 2F
convertFormat(nk.merged, from="seurat", to="anndata", outFile='nk_input.h5ad')
# Python in Spyder in Anaconda
import scanpy as sc
adata =sc.read_h5ad("nk_input.h5ad")
markers={ 'CYTOTOXICITY': ['GZMA', 'GZMB', 'GZMK', 'PRF1'], 'INHIBITORY RECEPTORS': ['LAIR1', 'KLRC1', 'TNFSF10', 'FASLG'], 'ACTIVATING RECEPTORS': ['FCGR3A', 'FCER1G', 'KLRC2', 'KLRK1'], 'CYTOKINES AND CHEMOKINES': ['IFNG', 'LTB', 'CCL3', 'XCL2'], 'CHEMOKINE AND CYTOKINE RECEPTORS': ['IL21R', 'CXCR3', 'CXCR4', 'CCR7'],  'ADHESION MOLECULES': ['ITGAL', 'CD58', 'SELL', 'CD96'], }
sc.pl.stacked_violin(adata, markers, 'cluster', dendrogram=False, save='.pdf')

# FIGURE 2G
convertFormat(nk.merged, from="seurat", to="anndata", outFile='nk_input.h5ad')
# Python in Spyder in Anaconda
import scanpy as sc
adata =sc.read_h5ad("nk_input.h5ad")
markers={ 'IL2': ['JAK3', 'STAT5A', 'IL2RA', 'SYK', 'PTK2B'], 'IL12': ['RALA', 'IL12RB2', 'IFNG', 'TCP1', 'SOD2'], 'IL15': ['STAT5', 'GRB2', 'IL2RB', 'IL2RG', 'IL15RA'], 'IL18': ['IL18R1', 'IL18RAP', 'IL18'], 'IL21': ['IL21R', 'STAT1', 'IL21'],  'STING': ['CGAS', 'DDX41', 'XRCC5', 'PRKDC', 'TMEM173'], }
sc.pl.dotplot(adata, markers, 'number', dendrogram=False, cmap='seismic', standard_scale='var', save='.pdf')

# Supplementary Figure 1B
oncotice.markers <- FindAllMarkers(oncotice, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
oncotice.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(oncotice, features = top10$gene, group.colors = col_vector, size = 2, group.bar.height = 0.01) + NoLegend() + scale_fill_gradientn(colors = c('white',"white", "#D33682"))

# Supplementary Figure 2C
nk.merged.markers <- FindAllMarkers(nk.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nk.merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(nk.merged, features = top10$gene, group.colors = c('#457b9d', '#ffb942', '#e63946', 'gray40', 'gray60', 'gray80'), size = 2, group.bar.height = 0.01) + NoLegend() + scale_fill_gradientn(colors = c('white',"white", "#D33682"))

# Supplementary Figure 2D
gene.sets <- getGeneSets(library = "C2", gene.sets = 'REACTOME_INTERLEUKIN_2_SIGNALING')
ES <- enrichIt(obj = nk.merged, gene.sets = gene.sets)
nk.merged <- AddMetaData(nk.merged, ES)
ES2 <- data.frame(nk.merged[[]], Idents(nk.merged))
colnames(ES2)[ncol(ES2)] <- "number"
ES2 <- ES2[,27:54]
ridgeEnrichment(ES2, gene.set = "REACTOME_INTERLEUKIN_2_SIGNALING", group = "number", add.rug = TRUE)
