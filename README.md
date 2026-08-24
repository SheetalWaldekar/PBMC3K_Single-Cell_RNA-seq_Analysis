## PBMC3K Single-Cell RNA-seq Analysis with R and Seurat

## Project Overview

This project performs an end-to-end single-cell RNA sequencing (scRNA-seq) analysis of the publicly available PBMC3K dataset from 10x Genomics.

The dataset contains approximately 2,700 human peripheral blood mononuclear cells (PBMCs) and was generated using 10x Genomics technology.

## Objective

The goal is to process the expression matrix, identify transcriptionally distinct cell populations, visualize the cells, identify cluster-specific marker genes, and annotate the clusters using known immune-cell markers.

## Biological Question

What major immune cell populations can be identified from the PBMC3K single-cell RNA-seq dataset based on their gene-expression profiles?

## Dataset

Dataset: 10x Genomics PBMC3K
Organism: Homo Sapians
Initial Cells: 2700
Initial genes: 32738
Technology: 10x Genomics single-cell RNA-seq
Genome annotation in the supplied matrix: hg19
Download source: 10x Genomics PBMC3K dataset

The raw 10x matrix is intentionally not required in this GitHub repository. It can be downloaded from the original source using the shell script.

## Tools and Technologies

Linux / Ubuntu
R 
Seurat
ggplot2
dplyr
patchwork
PCA 
UMAP
Marker-gene analysis

## Analysis Workflow

10x Genomics expression matrix
           ↓
     Data loading
           ↓
     Create Seurat Object
           ↓
     Quality control
           ↓
       Cell/gene filtering
           ↓
       Normalization
           ↓
 Highly variable gene selection
           ↓
       Scale Data
           ↓
          PCA
           ↓
  Nearest-neighbor graph
           ↓
      clustering
           ↓
         UMAP 
           ↓         
   Marker gene analysis
           ↓
   Cell-type annotation
           ↓
   Biological Interpritation

## Result

After quality-control filtering, 2640 cells remained from the original PBMC3K dataset of approximately 2700 cells and 13714 genes remained from aproximately 32738.
PCA-based dimensionality reduction followed by graph-based clustering identified 9 transcriptionally distinct cell clusters.
The clusters contained 684, 481, 476, 344, 291, 162, 155, 32, and 13 cells, respectively.
Cell-type annotation identified nine distinct populations: naive CD4 T cells, memory CD4 T cells, CD8 T cells, CD14+ classical monocytes, FCGR3A+ monocytes, B cells, NK cells, dendritic cells, and platelets. These populations represent major lymphoid, myeloid, antigen-presenting, and platelet lineages present in peripheral blood.

Cluster	Cell type:
Cluster 0         Naive CD4 T
Cluster 1	CD14+ Monocytes
Cluster 2	Memory CD4 T
Cluster 3	B cells
Cluster 4	CD8 T
Cluster 5	FCGR3A+ Monocytes
Cluster 6	NK cells
Cluster 7	Dendritic cells
Cluster 8	Platelets

## Conclusion

This project demonstrates a complete scRNA-seq analysis workflow from a raw 10x Genomics expression matrix to biologically interpretable cell populations.
The analysis demonstrates practical skills in:
Single-cell RNA-seq preprocessing through Quality control, Normalization, Dimensionality reduction, Clustering, Marker gene analysis, Cell-type annotation, Biological interpretation, Reproducible computational analysis

## Major Cell-Type Markers

T cells: CD3D, CD3E
B cells: MS4A1, CD79A, CD79B, CD37
NK cells: NKG7, GNLY, PRF1, GZMA, GZMM, CCL5
Monocytes: LYZ, LST1, TYROBP, AIF1, FCN1, S100A8, S100A9
Dendritic cells: FCER1A, CST3, CD74, HLA-DRA, HLA-DPA1, HLA-DPB1
Platelet/Megakaryocyte-like cells: PF4, PPBP, RGS18, SDPR

## Visualization and Interpretation

**QC Violin Plot:** Shows the distribution of genes detected per cell, total RNA counts per cell, and mitochondrial RNA percentage.
**PCA Plot:** Shows the major sources of transcriptional variation captured by principal components.
**Elbow Plot:** Helps determine the number of principal components used for downstream analysis.
**UMAP Cluster Plot:** Displays computationally identified cell populations. Each point represents a cell.
**Marker DotPlot:** Compares marker-gene expression across clusters. Dot size represents the percentage of cells expressing the gene, while dot intensity represents average expression.
**FeaturePlot:** Shows where specific marker genes are expressed on the UMAP.

## Figures

**01_QC_violin.png**: Distribution of genes, counts and mitochondrial percentage
**02_QC_scatter.png:** Relationships between QC metrics
**03_variable_features.png:** Highly variable genes
**04_PCA.png:** Principal component representation
**05_ElbowPlot.png:** Evaluation of principal components
**06_UMAP.png:** UMAP colored by computational clusters
**07_marker_DotPlot.png:** Canonical marker expression across clusters
**08_FeaturePlot.png:** Marker-gene expression on UMAP
**09_Final_annotated_UMAP.png:** Final biological cell-type annotation

## Author

Sheetal

