library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(hdf5r)

# Load the PBMC dataset
pbmc.data <- Read10X_h5("data/10k_PBMC.h5")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)

# Add mitochondrial content as a percentage
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells in the control group
head(pbmc@meta.data, 5)

#Normalize data
pbmc <- NormalizeData(pbmc)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2


# Here, we filter away cells that have unique feature counts(genes) over 5,000 or less than 200. 
# We also filter away cells that have > 15% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

#We can visualize QC metrics again after filtering cells
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2

# Depend on the data we analyze, we can use different cutoff for 
# nFeature_RNA and percent.mt. For example, we can filter away cells 
# that have unique feature counts(genes) over 5,000 or less than 300, or cells 
# that have > 10% mitochondrial counts, and see how QC metrics looks like.
temp <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
VlnPlot(temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.001)

plot1 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2

rm(temp) #removing unwanted cells from the dataset

# By default, we employ a global-scaling normalization method “LogNormalize” that 
# normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
#Normalized values are stored in pbmc[["RNA"]]@data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# set seed and put two plots in one figure
set.seed(123)
par(mfrow=c(1,2))
# original expression distribution
raw_geneExp = as.vector(pbmc[['RNA']]$counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp)
# expression distribution after normalization
logNorm_geneExp = as.vector(pbmc[['RNA']]$data) %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
hist(logNorm_geneExp) 

# same as the above command, skip here
pbmc <- NormalizeData(pbmc)

#We next calculate a subset of features 
#that exhibit high cell-to-cell variation in the dataset 
#(i.e, they are highly expressed in some cells, and lowly expressed in others
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")
plot1 + plot2

#Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing 
#step prior to dimensional reduction techniques like PCA. The ScaleData() function:
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose = FALSE)
#When using the above command, we use all genes to scale data. Scaling is an essential step 
#in the Seurat workflow, but only on genes that will be used as input to PCA. Therefore, 
#the default in ScaleData() is only to perform scaling on the previously identified variable 
#features (2,000 by default). And it will make this step faster.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
#Next we perform PCA on the scaled data. By default, only the previously determined variable features
#are used as input, but can be defined using features argument if you wish to choose a different subset.
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)

pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
#K-nearest neighbor (KNN) graph,with edges drawn between cells with similar feature expression patterns
# dims:20 = first 20 principle components
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
#To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default)
#or SLM (Blondel et al. 2008), to iteratively group cells together, with the goal of optimizing the standard modularity function. 
Idents(pbmc) #  The clusters can be found

head(Idents(pbmc), 5)
# Look at cluster IDs of the first 5 cells

#Seurat offers several non-linear dimensional reduction techniques, 
#such as tSNE and UMAP, to visualize and explore these datasets.
pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)
#Then we can get the UMAP plot of the single cell clustering results.
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:20, verbose = FALSE)
DimPlot(pbmc, reduction = "tsne")

DimPlot(pbmc, reduction = "tsne", label = TRUE)

plot <- DimPlot(object = pbmc)
LabelClusters(plot = plot, id = 'ident')

saveRDS(pbmc, file = "C:/Users/tubak/OneDrive/Belgeler/R_workspace/pbmc_processed.rds")



