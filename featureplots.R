# Load necessary libraries
library(dplyr)
library(Seurat)
library(ggplot2)

# Read the data
data <- read.csv("/N/scratch/syennapu/sigclust_siggenes.csv", header = TRUE, sep = '\t')

#View(data)
# Define a scoring function to prioritize markers
# Adjust the weights for each metric as needed
score_marker <- function(df) {
  df %>%
    mutate(score = log2(rank(-p_val_adj)) + log2(rank(abs(avg_log2FC)) )) %>%
    arrange(desc(score))
}

# Apply the scoring function to each cluster and select the top 5 markers
top_genes <- data %>%
  group_by(cluster) %>%
  do(score_marker(.)) %>%
  slice_max(order_by = score, n = 5) %>%
  ungroup()

# Print the top genes for each cluster
print(top_genes)
# Ensure the top genes are unique across clusters
top_genes <- top_genes %>%
  distinct(gene, .keep_all = TRUE)

# Print the top genes for each cluster
print(top_genes)
View(top_genes)
# Load your single-cell data object (assuming it's in a Seurat object format)
data1 <- readRDS("/N/scratch/syennapu/seurat_object.rds")
data1 <- FindNeighbors(data1, dims = 1:15)
data1 <- FindClusters(data1, resolution = 0.5)
head(Idents(data1), 5)
data1 <- RunUMAP(data1, dims = 1:15)
# Check if the genes are present in the Seurat object
available_genes <- rownames(data1)

missing_genes <- top_genes$gene[!top_genes$gene %in% available_genes]
if (length(missing_genes) > 0) {
  warning(paste("The following genes are not found in the Seurat object:", paste(missing_genes, collapse = ", ")))
} else {
  for (gene in unique(top_genes$gene)) {
    p <- FeaturePlot(data1, features = gene, raster= F , pt.size = 0.5) + 
      ggtitle(paste("Feature Plot for Gene:", gene)) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
}