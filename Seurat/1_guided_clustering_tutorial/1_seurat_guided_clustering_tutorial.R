## Tutorial from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html


## A. INSTALLATION
## B. READ-IN
## C. QC CONTROL/METRICS
## D. NORMALIZING DATA
## E. IDENTIFICATION OF HIGHLY VARIABLE FEATURES
## F. SCALING THE DATA
## G. LINEAR DIMENSIONAL REDUCTION
## H. DETERMINE DIMENSIONALITY OF THE DATASET
## I. CLUSTER THE CELLS
## J. RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)
## K. CLUSTER BIOMARKERS
## L. ASSIGNING CELL TYPE TO CLUSTERS


## A. INSTALLATION - using RStudio (https://posit.co/downloads/)
## A. Installation - using RTools (https://cran.rstudio.com/bin/windows/Rtools/)

## A. Enter commands in R (or R studio, if installed)
# install.packages('Seurat')
# install.packages('reticulate')
library(Seurat)


## A. Dataset being analyzed are PBMCs from 10x genomics 
## A. (2700 single cells sequenced on Illumina NextSeq 500)
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)


## B. READ-IN
## B. Read10X() fucntion reads in output of cellranger pipeline, returns unique 
##    molecular identified (UMI) count matrix. 

## B. Set working directory
setwd("D:/Dropbox/scRNAseq/scRNAseq_tutorials/Seurat/1_guided_clustering_tutorial/data")

## B. Read in data and store in pbmc.data variable
pbmc.data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")

## B. Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

## B. Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

## B. Dot values in matrix represent 0s - no molecules detected. Seurat 
##    represents these as a . value to save space (sparse-matrix)

## B. Dense size (with 0's)
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

## B. Sparse size (with .'s)
sparse.size <- object.size(pbmc.data)
sparse.size

## B. Ratio of dense to sparse
dense.size/sparse.size


## C. QC CONTROL/METRICS
## 1. Number of unique genes detected in each cell
##    1.1 Low-quality cells or empty droplets have LOW gene count
##    1.2 Cell doublets or multiplets may have HIGH gene count
## 2. Total number of molecules within cell correlates with unique genes
## 3. Percentage of reads that map to mitochondrial genome
##    3.1 Low-quality/dying cells exhibit extensive mitochondrial contamination
##    3.2 Mitochondrial QC metrics calculatated with PercentageFeatureSet() 
##        function calculates percentage of counts originating from a set of 
##        features
##    3.3 Use the set of all genes starting with MT - as a set of mitochondrial 
##        genes

## C. The [[ operator can add columns to object metadata. This is a great place 
##    to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## C. QC metrics - number of unique genes and total molecules are automatically 
##    calculated during CreateSeuratObject() and stored in object metadata

## C. Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

##C. Visualize QC metrics and use these to filter cells
##   - Filter cells that have unique feature counts over 2500 (doublets/multi) 
##     or less than 200 (low-quality cells)
##   - Filter cells that have >5% mitochondrial counts

## C. Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## C. FeatureScatter used to visualize feature-feature relationships but can be
##    used for anything calculated by object - columns: metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## C. Select PBMC subset that meets criteria above
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## D. NORMALIZING DATA
## D. After removing unwanted cells, need to normalize the data. Employ a 
##    global-scaling normalization method: "LogNormalize" that normalizes the
##    feature expression measurements for each cell by the total expression,
##    multiplies this by a scale factor (10,000 by default) and log-transforms

## D. Normalized values are stored in pbmc[["RNA]]#data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
## N.B. same behaviour achieved by pbmc <- NormalizeData(pbmc) (default values
## for function showed above)


## E. IDENTIFICATION OF HIGHLY VARIABLE FEATURES
## E. Calculate a subset of features that exhibit high cell-to-cell variation
##    (i.e. highly expressed in some cells and lowly expressed in others).
##    Focusing on these genes in downstream analysis helps to highlight biological
##    signal in single-cell datasets

## E. Models mean-variance relationship inherent in single-cell data and is 
##    implemented in FindVariableFeatures() function. By default, they return
##    2,000 features per dataset - used in downstream analysis 
##    (e.g. Principal Component Analysis - PCA)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## E. Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

## E. plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


## F. SCALING THE DATA
## F. Apply a linear transformation ('scaling') that is a standard pre-process
##    step prior to dimensional reduction techniques like PCA
## F. ScaleData() function:
##    1. Shifts expression of each gene so that mean expression across cells is 0
##    2. Scales expression of each gene so that variance across cells is 1
##       - Step gives equal weight in downstream analyses so that highly-
##         expressed genes do not dominate
##    3. Results are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
## N.B. Default in ScaleData() is to perform scaling only on variable features
## that will be used as input to PCA. Can leave out all.genes argument in
## prior function call if takes a long time. Seurat heatmaps require genes in
## heatmap to be scaled to make sure highly-expressed genes dont dominate the 
## heatmap - all gense are being scaled in this tutorial

## N.B. SCTransform() in Seurat v3 can regress out unwanted variation which
## might be associated with cell cycle stage, or mitochondrial contamination
## (https://satijalab.org/seurat/articles/sctransform_vignette.html)


## G. LINEAR DIMENSIONAL REDUCTION
## G. Now perform principal component analysis on the scaled data. By default, 
##    only previously determined variable features are used as input, but can 
##    be defined using features argument if you use a different subset
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

## G. Seurat provides several useful ways of visualizing both cells and features
##    that define the PCA, including VizDimReduction(), DimPlot(), and 
##    DimHeatmap()

## G. Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

## G. DimHeatmap() allows for easy exploration of the primary sources of
##    heterogeneity in a dataset, and can be useful when trying to decide which
##    PCs to include for further downstream analyses.

## G. Both cells and features are ordered according to their PCA scores. Setting
##    cells to a number plots the 'extreme' cells on both ends of the spectrum,
##    which dramtically speeds plotting for large datasets.
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


## H. DETERMINE DIMENSIONALITY OF THE DATASET
## H. To overcome extensive technical noise in any single feature for data,
##    Seurat clusters cells based on their PCA scores, with each PC representing
##    a metafeature that combines information across a correlated feature set

## H. Top principal components represent a robust compression of the dataset

## H. How to determine how many components to include? JackStraw plot - wherein
##    randomly permute a subset of the data (1% by default) and rerun PCA,
##    constructing a 'null distribution' of feature scores, and repeat procedure

## H. Significant PCs are those ID'd as having strong enrichment of low p-values

## H. NOTE: This process can take a long time for big datasets, comment out for 
##    expediency. More approximate techniques such as those implemented in 
##    ElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

## H. Can visualize distribution of p-values for each PC with a uniform
##    distribution (dashed line). 'Significant' PCs will show a strong 
##    enrichment of features with low p-values (solid curve above dashed line)
JackStrawPlot(pbmc, dims = 1:15)
## N.B. Sharp drop off in significance after first 10-12 PCs

## H. An alternative heuristic method generates an "Elbow plot': ranking of PCs
##    based on percentage of variance explained by each one (ElbowPlot() function)
ElbowPlot(pbmc)
## N.B. Observe an 'elbow' around PC9-10, suggesting that majority of true 
##      signal captured in first 10 PCs.

## H. Identifying true dimensionality with rarer populations can be hard. 3 ways:
##    1. Exploring PCs to determine sources of heterogeneity, use with GSEA
##    2. Implements statistical test based on random null model
##    3. Heuristic that can be commonly used and calculated instantly

## H. Chose 10 items here - can try repeat downstream analyses with different
##    numbers of PCs (10,15 or even 50!). Err on higher side when possible


## I. CLUSTER THE CELLS
## I. Graph-based clustering approach. Cells are embedded in a graph structure
##    for example a K-nearest neighbor graph (KNN) with edges drawn between cells
##    with similar feature expression patterns, and then attempt to partition into
##    highly interconnected 'quasi-cliques' or 'communities'.

## I. Determine distance in space using FindNeighbors() function. Use the
##    dimensionality of the dataset determined above (first 10 PCs)

## I. Iteratively group cells together, with goal of optimizing standard 
##    modularity function. FindClusters() function implements this procedure and
##    contains a resolution paramater that sets the granularity of the
##    downstream clustering, with increased values leading to greater number of
##    clusters. Setting this parameter between 0.4-1.2 returns good results for 
##    single-cell datasets of around 3K cells. Optimal resolution increases for 
##    larger datasets. Clusters can then be found using Idents() function.
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

## I. Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


## J. RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)
## J. tSNE and UMAP visualizes and explores these datasets. Places similar cells
##    together in low-dimensional space. Use same PCs as input to clustering
##    analysis (10)

## J. If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
## 'umap-learn')
# reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

## J. Note that you can set `label = TRUE` or use the LabelClusters function to 
##    help label individual clusters
DimPlot(pbmc, reduction = "umap")

## J. Can save the object at this point so it can be easily loaded back in 
##    without having to rerun the computationally intensitve steps above, or
##    easily shared with collaborators
dir.create("/output/")
saveRDS(pbmc, file = "D:/Dropbox/scRNAseq/scRNAseq_tutorials/Seurat/1_guided_clustering_tutorial/output/pbmc_tutorial.rds")

## N.B. To load saved RDS file back into R: 
##      pbmc <- readRDS(file = "D:/Dropbox/scRNAseq/scRNAseq_tutorials/Seurat/1_guided_clustering_tutorial/output/pbmc_tutorial.rds")


## K. CLUSTER BIOMARKERS
## K. Seurat helps to find markers that define clusters via differential exp.
##    By default, it identified positive and negative markers of a single cluster
##    (specificed in ident.1), compared to all other cells

## K. FindAllMarkers() automates this process for all clusters, but can test
##    groups of clusters vs each other or against all cells

## K. min.pct argument requires a feature to be detected at a minimum percentage
##    in either of the two groups of cells and the thresh.test argument requires
##    a feature to be differentially expressed (on average) by some amount between
##    the two groups. Can set these both to 0, but with big increase in time - will
##    set up a large number of features that are unlikely to be highly discriminatory

## K. Another option to speed up these computations, max.cells.per.ident can be set
##    This will downsample each identity class to have no more cells than whatever
##    this is set to. While this will result in a loss in power, speed increase is
##    significant and most highly differentially expressed features will rise to top.


## K. Find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

## K. Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

## K. Find markers for every cluster compared to all remaining cells, report only the positive
##    ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

## K. Seurate has several tests for differential expression that can be set with
##    the test.use parameter. ROC test returns the classification power for any
##    individual marker (ranging from 0-random to 1-perfect)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

## K. Several tools to visualize marker expression. VlnPlot() shows expression
##    probability distributions across clusters and FeaturePlot() visualizes
##    feature expression on a tSNE or PCA plot are most commonly used visualizations

## K. Can also try RidgePlot(), CellScatter(), DotPlot()
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

## K. You can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

## K. DoHeatmap() generates an expression heatmap for given cells and features. 
##    In this case, we are plotting the top 20 markers (or all markers if less 
##    than 20) for each cluster.
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


## L. ASSIGNING CELL TYPE TO CLUSTERS
##    Cluster ID	Markers	     Cell Type
##    0	          IL7R, CCR7	 Naive CD4+ T
##    1	          CD14, LYZ	   CD14+ Mono
##    2	          IL7R,S100A4  Memory CD4+
##    3	          MS4A1	       B
##    4	          CD8A	CD8+   T
##    5	          FCGR3A,MS4A7 FCGR3A+ Mono
##    6	          GNLY, NKG7	 NK
##    7	          FCER1A,CST3	 DC
##    8	          PPBP	       Platelet
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "D:/Dropbox/scRNAseq/scRNAseq_tutorials/Seurat/1_guided_clustering_tutorial/output/pbmc3k_final.rds")
