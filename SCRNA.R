install.packages("Seurat")
install.packages(c("ggplot2", "dplyr", "patchwork"))
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

## Load data
pbmc3k_data <- Read10X(
  data.dir= "//wsl.localhost/Ubuntu-24.04/home/shital/my_project/scrna_project/filtered_gene_bc_matrices/hg19"
)
dim (pbmc3k_data)

## create seurat Object
pbmc <- CreateSeuratObject(counts = pbmc3k_data, project = 'PBMC3K', min.cells = 3, min.features = 200)
pbmc

# QC
# mitochondrial % calculate
pbmc[['mt_percent']] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize QC Matrix
VlnPlot(pbmc, features = c('nCount_RNA', 'nFeature_RNA', 'mt_percent' ))

#Counts vs Mitochondrial %
plot1<- FeatureScatter(pbmc, feature1='nCount_RNA', feature2='mt_percent')

# Counts vs genes
plot2<- FeatureScatter(pbmc, feature1='nCount_RNA', feature2='nFeature_RNA')
plot1 + plot2

# Filter Low quality cells
pbmc <- subset(pbmc, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & mt_percent < 5)
pbmc

#Normalization and feature selection
pbmc <- NormalizeData(pbmc)
pbmc

# find highly variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)
plot1 <- VariableFeaturePlot(pbmc)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale Data (prepare Data for PCA plot)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc

#PCA plot
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
print(pbmc[["pca"]], dims = 1:5)

#visualize top genes driving each pc
VizDimLoadings(pbmc, dims=1:2, reduction='pca')

DimPlot(
  pbmc,
  reduction = "pca"
)

# Elbow plot to Decide how many PCs to use
ElbowPlot(pbmc)

#find Nearest neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc

#Clustering
pbmc <- FindClusters(pbmc,resolution = 0.5)
print(paste("Number of clusters:", length(levels(pbmc))))

table(Idents(pbmc))

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

#first major biological figure
DimPlot(pbmc, reduction = 'umap', label=TRUE)

# Find markers for every cluster vs all others
pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(pbmc_markers)

# Top markers per cluster
pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by=avg_log2FC)

# Heatmap of top 10 markers per cluster
top10 <- pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
top10
DoHeatmap(pbmc, features=top10$gene) + NoLegend()

# DotPlot Visualize marker genes 
markers_to_plot <- c("CD3D", "CD3E", "MS4A1", "CD79A", "NKG7", "GNLY", "LYZ", "LST1", "FCER1A", "CST3", "PF4", "PPBP")

DotPlot(pbmc, features = markers_to_plot) + RotatedAxis()

# Visualize known markers
FeaturePlot(pbmc, features=c('CD3D','CD4','CD8A','MS4A1', 'GNLY','CD14','LYZ','FCGR3A'), ncol=3)

# Assign cell type labels (adjust cluster numbers to match yours)
new_cluster_ids <- c('Naive CD4 T', 'CD14+ Mono', 'Memory CD4 T', 'B', 'CD8 T', 'FCGR3A+ Mono', 'NK', 'DC', 'Platelet')
new_cluster_ids
names(new_cluster_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_cluster_ids)

# Final annotated UMAP
DimPlot(pbmc, reduction='umap', label=TRUE, pt.size=0.5) + NoLegend()

#Compare CD14 Monocytes vs CD16 Monocytes
mono.markers <- FindMarkers(pbmc, ident.1='CD14+ Mono', ident.2='FCGR3A+ Mono', min.pct=0.25)
head(mono.markers, 10) 

saveRDS(pbmc, file='pbmc3k_analyzed.rds')
pbmc <- readRDS('pbmc3k_analyzed.rds')

write.csv(Idents(pbmc), file='cell_type_assignments.csv')










































