install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")  # Human gene annotation package
install.packages("dplyr")   # For data manipulation
install.packages("ggplot2")  # For visualization
install.packages("pheatmap") # For heatmap visualization
install.packages("VennDiagram")
install.packages("openxlsx") # For saving to Excel
# Load required libraries
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble)
library(dplyr)
library(ggplot2)
library(pheatmap) # For heatmap visualization
library(VennDiagram)
library(openxlsx)   # For saving to Excel
# Set your working directory to the folder where your files are located
setwd("C:\\Users\\hp\\Downloads\\limma")  
getwd()

# Load the data
table1 <- read.table("GSE164332.top.table (1).tsv", sep = "\t", header = TRUE)
table2 <- read.table("GSE122063.top.table (3).tsv", sep = "\t", header = TRUE)

# Standardize column names for merging
table1 <- table1 %>% rename(Gene = Symbol, log2FC_1 = log2FoldChange, p_value_1 = pvalue, padj_1 = padj)
table2 <- table2 %>% rename(Gene = ID, log2FC_2 = logFC, p_value_2 = P.Value, padj_2 = adj.P.Val)

#After renaming columns in table2
table2$GeneSymbol <- mapIds(org.Hs.eg.db, keys = table2$GB_ACC, column = "SYMBOL", keytype = "ACCNUM", multiVals = "first")

# Merge the two tables by Gene
merged_table <- merge(table1, table2, by.x = "Gene", by.y = "GeneSymbol")
table2$GeneSymbol <- mapIds(org.Hs.eg.db, 
                            keys = table2$GB_ACC, 
                            column = "SYMBOL", 
                            keytype = "ACCNUM", 
                            multiVals = "first")

# Find shared genes (intersection of Gene symbols)
shared_genes <- intersect(table1$Gene, table2$GeneSymbol)

# Subset both tables for only shared genes
table1_shared <- table1 %>% filter(Gene %in% shared_genes)
table2_shared <- table2 %>% filter(GeneSymbol %in% shared_genes)

# Merge shared genes into a single dataframe
shared_genes_df <- merge(table1_shared, table2_shared, by.x = "Gene", by.y = "GeneSymbol")

# Optional: Add up/down regulation column based on log2FC in table1
shared_genes_df <- shared_genes_df %>%
  mutate(Regulation = case_when(
    log2FC_1 > 1.1 ~ "Upregulated",
    log2FC_1 < -1.1 ~ "Downregulated",
    TRUE ~ "Neutral"))   # closing the mutate() function

#Saving the tables of mutate function
write.csv(shared_genes_df, "shared_genes.csv", row.names = FALSE)
write.csv(shared_genes, "shared_genes 1.csv", row.names = FALSE)

# Filter based on log2 fold change and create additional columns for up/down regulation
# Check what the Gene column looks like in each
head(table1$Gene)
head(table2$Gene)


# Merge the two tables by Gene (now using gene symbols)
merged_table <- merge(table1, table2, by = "Gene")

# Check the merged table
head(merged_table)
nrow(merged_table)

# Filter based on log2 fold change and creating additional column for up/down regulation
filtered_DEGs <- merged_table %>%
  filter(padj_1 < 0.05 & padj_2 < 0.05 & 
           ((log2FC_1 > 1 & log2FC_2 > 1) | (log2FC_1 < -1 & log2FC_2 < -1))) %>%
  mutate(regulation = ifelse(log2FC_1 > 1 & log2FC_2 > 1, "Upregulated", "Downregulated"))

# Check filtered results
nrow(filtered_DEGs)

# Save the filtered results to file
write.csv(filtered_DEGs, "filtered_DEGs.csv", row.names = FALSE)

# Print first few rows
print(head(filtered_DEGs))
# Volcano Plot
ggplot(merged_table, aes(x = log2FC_1, y = -log10(padj_1))) +
  geom_point(aes(color = ifelse(log2FC_1 > 1.1, "Upregulated",
                                ifelse(log2FC_1 < -1.1, "Downregulated", "Neutral")))) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Neutral" = "gray")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)") +
  theme_minimal()

# Generate the heatmap with proper color mapping and gene selection
# Prepare data for heatmap: Select genes and log2 fold change
heatmap_data <- filtered_genes %>%
  dplyr::select(Gene, log2FC_1, log2FC_2)
rownames(heatmap_data) <- heatmap_data$Gene
heatmap_data <- heatmap_data %>%
  dplyr::select(-Gene) %>%  # Remove the Gene column
  as.matrix()


# Assign colors based on fold change for heatmap
heatmap_colors <- ifelse(heatmap_data[, 1] > 1.1, "red", 
                         ifelse(heatmap_data[, 1] < -1.1, "blue", "gray"))

# Define a custom color palette for the heatmap to match the volcano plot colors
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Load the tibble package again to use column_to_rownames
library(tibble)

# Prepare data for heatmap: Select genes and log2 fold change
# Select columns and set rownames manually
selected_data <- filtered_genes %>%
  dplyr::select(Gene, log2FC_1, log2FC_2)

rownames(selected_data) <- selected_data$Gene  # Set rownames
selected_data <- selected_data[, -1]  # Remove the "Gene" column
heatmap_data <- as.matrix(selected_data)  # Convert to matrix

# Assign colors based on fold change for heatmap
heatmap_colors <- ifelse(heatmap_data[, 1] > 1.1, "red", 
                         ifelse(heatmap_data[, 1] < -1.1, "blue", "gray"))

# Define a custom color palette for the heatmap to match the volcano plot colors
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Generate the heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         scale = "none",  # No scaling to show raw values, similar to volcano plot
         color = color_palette,  # Use the custom color palette
         main = "Heatmap of Log2 Fold Changes")

# Save the heatmap as a file
ggsave("heatmap_plot.png", width = 6, height = 6)

# Sort the res.df by log2FoldChange in descending order
up.df.sorted <- filtered_genes[order(-filtered_genes$log2FC_1), ]
down.df.sorted <- filtered_genes[order(-filtered_genes$log2FC_2), ]

# Save the sorted data frame as a text file
write.table(up.df.sorted, file = "up_res_df.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(down.df.sorted, file = "down_res_df.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# venn diagram preparation
genes_ds1 <- table1$Gene[table1$padj_1 < 0.05 & abs(table1$log2FC_1) > 1.1]
genes_ds2 <- table2$GeneSymbol[table2$padj_2 < 0.05 & abs(table2$log2FC_2) > 1.1]

# Remove NA values
genes_ds1 <- genes_ds1[!is.na(genes_ds1)]
genes_ds2 <- genes_ds2[!is.na(genes_ds2)]

# Print summary
cat("Dataset 1:", length(genes_ds1), "genes\n")
cat("Dataset 2:", length(genes_ds2), "genes\n")
cat("Shared:", length(intersect(genes_ds1, genes_ds2)), "genes\n")

# create venn diagram 
grid.newpage()
venn.plot <- venn.diagram(
  x = list(GSE164332 = genes_ds1, GSE122063 = genes_ds2),
  filename = NULL,  # NULL means display instead of save
  fill = c("orange", "purple"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  main = "DEGs Comparison"
)
grid.draw(venn.plot)
# Create and save Venn Diagram
venn.diagram(
  x = list(GSE164332 = genes_ds1, GSE122063 = genes_ds2),
  filename = "venn_diagram.png",
  fill = c("orange", "purple"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  main = "DEGs Comparison"
)