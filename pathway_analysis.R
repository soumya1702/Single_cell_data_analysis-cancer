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

# Plot KEGG enrichment results using adjusted p-value
kegg_plot <- ggplot(top_kegg, aes(x = reorder(Description, p.adjust), y = p.adjust)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  scale_size_continuous(range = c(2, 7)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "KEGG Pathway Enrichment Analysis",
       x = "Pathway",
       y = "Adjusted P-value",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(kegg_plot)



