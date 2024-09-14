# Load necessary libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggpubr)
# Read the data
data_path <- "/N/scratch/syennapu/sigclust_siggenes.csv"
all_data <- read.csv(data_path, header = TRUE, sep = '\t')

# Print the first few rows to understand the structure
head(all_data)
# Extract significant genes based on p_val_adj
significant_genes <- all_data %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# Convert gene symbols to Entrez IDs
gene_list <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = gene_list$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

# Extract the top KEGG pathways based on adjusted p-value
top_kegg <- kegg_enrichment@result %>% top_n(20, wt = -p.adjust)

# Remove "Mus musculus (house mouse)" from the Description
top_kegg$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", top_kegg$Description)


# Split geneID into individual genes for each pathway
top_kegg_genes <- top_kegg %>%
  separate_rows(geneID, sep = "/") %>%
  group_by(Description) %>%
  summarize(Genes = paste(unique(geneID), collapse = ", "))

# Convert Entrez IDs back to gene symbols
entrez_ids <- unlist(strsplit(top_kegg_genes$Genes, ", "))
entrez_to_symbol <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)

# Create a mapping of Entrez IDs to gene symbols
id_to_symbol <- setNames(entrez_to_symbol$SYMBOL, entrez_to_symbol$ENTREZID)

# Replace Entrez IDs with gene symbols in the summarized data frame
top_kegg_genes <- top_kegg_genes %>%
  rowwise() %>%
  mutate(Genes = paste(sapply(strsplit(Genes, ", ")[[1]], function(x) id_to_symbol[[x]]), collapse = ", "))

# Print the pathways and their associated genes
print(top_kegg_genes)

# Save the pathways and associated genes to a CSV file
write.csv(top_kegg_genes, "/N/scratch/syennapu/kegg_pathways_genes.csv", row.names = FALSE)
