# Advanced-Programming-and-Data-Analysis-Project-
Discovering the common pathways between SARS COV 2 and Alzheimer's disease

##Datasets

*GSE164332: COVID-19 gene expression data*

*GSE122063: Alzheimer's disease gene expression data*

## Requirements

Install the following R packages:

rBiocManager::install("org.Hs.eg.db")

install.packages(c("dplyr", "ggplot2", "pheatmap", "VennDiagram", "openxlsx"))

## What the Script Does

Loads and merges two gene expression datasets

Identifies shared genes between COVID-19 and Alzheimer's disease

Filters DEGs based on:

Adjusted p-value < 0.05

|log2 fold change| > 1


### Generates visualizations:

Volcano plot showing gene expression changes

Heatmap of log2 fold changes

Venn diagram comparing DEGs between datasets


## Output Files

`shared_genes.csv` - All genes present in both datasets

`filtered_DEGs.csv` - Significantly differentially expressed shared genes

`up_res_df.txt` / `down_res_df.txt` - Sorted upregulated/downregulated genes

`heatmap_plot.png` - Heatmap visualization

`venn_diagram.png` - Venn diagram of DEG overlap

## Usage

*1-Set your working directory to the folder containing your data files*

*2-Ensure GSE164332.top.table (1).tsv and GSE122063.top.table (3).tsv are in the directory*

*3-Run the script*

## Analysis Criteria

Upregulated: log2FC > 1

Downregulated: log2FC < -1

Significance: adjusted p-value < 0.05

# Previous project ( Covid vs. GBS )

Before the current project there was a project for discovering common pathways between covid (GSE157103) and GBS (GSE31014) 
but by coding on R we discovered that there were no shared genes but we discovered this after revising the code and we had already generated plots like heatmap and volcano but based on the wrong data and wrong filtration and missing analysis criteria.
