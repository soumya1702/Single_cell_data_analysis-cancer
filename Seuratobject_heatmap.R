library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load the dataset
data <- Read10X(data.dir = "/N/scratch/syennapu/GSE234577_RAW/")

# Initialize the Seurat object with the raw (non-normalized data).
d1 <- CreateSeuratObject(counts = data, project = "scRNA", min.cells = 3, min.features = 200)
d1
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
d1[["percent.mt"]] <- PercentageFeatureSet(d1, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(d1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(d1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
d1 <- subset(d1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalisation
d1 <- NormalizeData(d1, normalization.method = "LogNormalize", scale.factor = 10000)

d1 <- NormalizeData(d1)

d1 <- FindVariableFeatures(d1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(d1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(d1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(d1)
d1 <- ScaleData(d1, features = all.genes)

d1 <- RunPCA(d1, features = VariableFeatures(object = d1))

# Examine and visualize PCA results a few different ways
print(d1[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(d1, dims = 1:2, reduction = "pca")

DimPlot(d1, reduction = "pca") + NoLegend()

DimHeatmap(d1, dims = 1, cells = 500, balanced = TRUE)

#deciding dimensions based on pcs and elbowplot
ElbowPlot(d1)
DimHeatmap(d1, dims = 1:15, cells = 500, balanced = TRUE)

#finding number of neighbors and clusters
d1 <- FindNeighbors(d1, dims = 1:7)
d1 <- FindClusters(d1, resolution = 0.5)
head(Idents(d1), 5)
d1 <- RunUMAP(d1, dims = 1:7)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(d1, reduction = "umap", group.by = "seurat_clusters")

#saving as rds file
saveRDS(d1, file = "/N/scratch/syennapu/seurat_object.rds")
#d1 <- PrepSCTFindMarkers(d1)

#find all markers
d1_markers <- FindAllMarkers(d1,
                                   logfc.threshold = 0.25,
                                   min.pct = 0.1,
                                   only.pos = TRUE,
                                   recorrect_umi = F)
View(d1_markers)
output_directory <- "/N/u/syennapu/Quartz/"
output_file <- "allmarkers_d1.csv"
output_path <- file.path(output_directory, output_file)

#ensure the output directory exists
if (!dir.exists(output_directory)){
  dir.create(output_directory, recursive = TRUE)
}

#save the markers data frame to a specified directory
write.csv(d1_markers, output_path, row.names = FALSE)
cat("markers data saved to", output_path)
# Extract top 5 markers per cluster
top5_d1 <- d1_markers %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_log2FC)

# Visualize top 5 markers per cluster
View(top5_d1)

output_directory <- "N/u/syennapu/Quartz/"
output_file <- "top5markers_d1.csv"
output_path <- file.path(output_directory, output_file)

#ensure the output directory exists
if (!dir.exists(output_directory)){
  dir.create(output_directory, recursive = TRUE)
}

#save the markers data frame to a specified directory
write.csv(top5_d1, output_path, row.names = FALSE)
cat("markers data saved to", output_path)

#plotting heatmap for top5 markers for each cluster
d1_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
options(repr.plot.width=25, repr.plot.height=25) 
DoHeatmap(subset(d1, downsample = 100),
          
          features = top5$gene)  + 
  theme(text = element_text(size = 10))

