install.packages("dplyr")   # For data manipulation
install.packages("ggplot2")  # For visualization
install.packages("pheatmap") # For heatmap visualization
install.packages("rlang")   # For plotting library.
install.packages("ggplot2") # For helper library that ggplot2 depends on.
install.packages("openxlsx") # For saving to Excel
# Load required libraries
library(dplyr)
library(ggplot2)
library(pheatmap)  # For heatmap visualization
library(rlang)
library(ggplot2)
library(openxlsx)   # For saving to Excel

# Set your working directory to the folder where your files are located
setwd("C:\\Users\\hp\\Downloads\\limma") # Replace with your folder path
getwd()

# Load the data
table1 <- read.table("GSE157103.top.table.tsv", sep = "\t", header = TRUE,quote = "")
table2 <- read.table("GSE31014.top.table.tsv", sep = "\t", header = TRUE,quote = "")


# Standardize column names for merging - using Symbol as the main gene identifier
table1 <- table1 %>% rename(Gene = Symbol, log2FC_1 = log2FoldChange, padj_1 = padj)
table2 <- table2 %>% rename(Gene = Gene.symbol, log2FC_2 = logFC, padj_2 = adj.P.Val)

# Check the tables
head(table1)
head(table2)

# Check the column names
colnames(table1)
colnames(table2)

# Check what the Gene column looks like in each
head(table1$Gene)
head(table2$Gene)

# See if there's any overlap at all
length(intersect(table1$Gene, table2$Gene))

# Merge the two tables by Gene (now using gene symbols)
merged_table <- merge(table1, table2, by = "Gene")

# Check the merged table
head(merged_table)
nrow(merged_table)

# Check the range of fold changes for table 1
summary(merged_table$log2FC_1)
min(merged_table$log2FC_1, na.rm = TRUE) # To decide the threshold
max(merged_table$log2FC_1, na.rm = TRUE) # To decide the threshold

summary(merged_table$log2FC_2)
min(merged_table$log2FC_2, na.rm = TRUE) # To decide the threshold
max(merged_table$log2FC_2, na.rm = TRUE) # To decide the threshold

summary(merged_table$padj_1)
min(merged_table$padj_1, na.rm = TRUE) # To decide the threshold
max(merged_table$padj_1, na.rm = TRUE)# To decide the threshold

summary(merged_table$padj_2)
min(merged_table$padj_2, na.rm = TRUE) # To decide the threshold
max(merged_table$padj_2, na.rm = TRUE)# To decide the threshold

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

# Load the tibble package to use column_to_rownames
library(tibble)

# Prepare data for heatmap: Select genes and log2 fold change
# Select columns and set rownames manually
selected_data <- filtered_genes %>%dplyr::select(Gene, log2FC_1, log2FC_2)

rownames(selected_data) <- selected_data$Gene  # Set rownames

selected_data <- selected_data[, -1]  # Remove the "Gene" column

heatmap_data <- as.matrix(selected_data)

# Assign colors based on fold change for heatmap
heatmap_colors <- ifelse(heatmap_data[, 1] > 1.1, "red", ifelse(heatmap_data[, 1] < -1.1, "blue", "gray"))

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

# Generate excel sheet of shared genes between both diseases
# Change gene names to uppercase to be the same and to avoid letter problems
genes1 <- toupper(table1$Gene)
genes2 <- toupper(table2$Gene)

# Shared genes
shared_genes <- intersect(genes1, genes2)
length(shared_genes)  # Number of genes
# Filter each table to get the specific common values between genes
shared_table1 <- table1 %>% filter(toupper(Gene) %in% shared_genes)
shared_table2 <- table2 %>% filter(toupper(Gene) %in% shared_genes)

# mergeing both tables into one based on genes
shared_merged <- merge(shared_table1, shared_table2, by = "Gene")
head(shared_merged)
nrow(shared_merged) 
write.csv(shared_merged, "shared_genes.csv", row.names = FALSE)