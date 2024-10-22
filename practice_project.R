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
                      

