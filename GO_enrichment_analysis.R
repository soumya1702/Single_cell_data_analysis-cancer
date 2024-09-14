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

# Create a dot plot for marker gene expression
dot_plot <- ggplot(all_data, aes(x = gene, y = as.factor(cluster))) +
  geom_point(aes(size = pct.1, color = avg_log2FC)) +
  theme_minimal() +
  labs(title = "Marker Gene Expression",
       x = "Marker Gene",
       y = "Cluster",
       size = "Fraction of Cells Expressing",
       color = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(dot_plot)

# Extract significant genes based on p_val_adj
significant_genes <- all_data %>%
  filter(p_val_adj < 0.05) %>%
  pull(gene)

# Convert gene symbols to Entrez IDs
gene_list <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform GO enrichment analysis for BP, CC, and MF
go_bp <- enrichGO(gene = gene_list$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
go_cc <- enrichGO(gene = gene_list$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
go_mf <- enrichGO(gene = gene_list$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# Extract top 5 GO terms for each category
top_go_bp <- go_bp@result %>% top_n(5, wt = -p.adjust)
top_go_cc <- go_cc@result %>% top_n(5, wt = -p.adjust)
top_go_mf <- go_mf@result %>% top_n(5, wt = -p.adjust)

# Combine the top GO terms
go_combined <- rbind(
  data.frame(Category = "BP", top_go_bp),
  data.frame(Category = "CC", top_go_cc),
  data.frame(Category = "MF", top_go_mf)
)

# Arrange the GO terms by Category and Description to ensure correct sequence
go_combined$Description <- factor(go_combined$Description, 
                                  levels = unique(go_combined$Description[order(go_combined$Category, go_combined$p.adjust)]))

# Plot combined GO enrichment results in a bar plot
go_combined_plot <- ggplot(go_combined, aes(x = Description, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 5 GO Enrichment Analysis",
       x = "GO Term",
       y = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(go_combined_plot)


