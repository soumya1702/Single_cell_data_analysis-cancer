d3  <- read.csv("/N/scratch/syennapu/pvals.csv",header = T, sep='\t')

d6 <- read.csv("/N/u/syennapu/Quartz/allmarkers_d1.csv")

# Identify significant clusters
# Define significance thresholds
significance_threshold <- 0.05
non_significance_threshold <- 0.05
significant_clusters <- rownames(d3)[d3$condition < significance_threshold & 
                                       d3$sex < non_significance_threshold & 
                                       d3$condition.sex < non_significance_threshold ]
                                       
# Remove the "Cluster " prefix from the significant clusters
significant_clusters <- gsub("Cluster ", "", significant_clusters)
print(significant_clusters)  # Check the adjusted cluster identifiers

# Subset the Seurat object for significant clusters
significant_data <- subset(d6, cluster %in% c('2','3','4'))
View(significant_data)


all <- subset(significant_data, avg_log2FC > 0.5 & p_val_adj < 0.05)
View(all)

write.table(all,"/N/scratch/syennapu/sigclust_siggenes.csv",row.names = T,col.names = T,sep='\t') 





