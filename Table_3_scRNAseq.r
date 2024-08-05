library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)

rm(list = ls())
memory.limit(size=56000)

# Load data
data <- LoadH5Seurat("521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat", assays = "data")

# Extract only "LD" data
data_ld <- subset(data, subset = sampletype == "LD")
data_ld$sampletype <- factor(data_ld$sampletype, levels = c("LD"))

# Set cluster identifiers
data_ld <- SetIdent(data_ld, value = data_ld@meta.data$subclass.l1)

# Save UMAP plot output in TIFF format
tiff("Umap_KPMP_LD.tiff", height=5, width=7, units="in", res=300)
DimPlot(data_ld, reduction = "umap", label = TRUE, repel = TRUE, order = c("Immune","Interstitial","IC","PC","CNT","DCT","TAL","ATL/TAL","DTL","PT","POD","PEC","EC"))
dev.off()

tiff("Umap_KPMP_LD_nolegend.tiff", height=5, width=5.2, units="in", res=300)
DimPlot(data_ld, reduction = "umap", label = FALSE, repel = TRUE) + NoLegend()
dev.off()

# Save UMAP plot output in TIFF format with legend and other elements removed
tiff("Umap_KPMP_LD_with_legend_dots.tiff", height=5, width=7, units="in", res=300)
DimPlot(data_ld, reduction = "umap", label = FALSE, repel = TRUE) +
  theme(
    legend.text = element_blank(),  
    legend.title = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(),  
    axis.text.y = element_blank()   
  )
dev.off()

# Save FeaturePlot output in TIFF format

tiff("feature_AKT1_data_ld.tiff", height=7, width=9, units="in", res=300)
FeaturePlot(data_ld, features = c("AKT1"))
dev.off()

# Save FeaturePlot output in TIFF format (without labels)
tiff("feature_AKT1_data_ld.tiff", height=7, width=9, units="in", res=300)
FeaturePlot(data_ld, features = c("AKT1"), label = FALSE) +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    plot.title = element_blank()  
  )
dev.off()

tiff("feature_SGK1_data_ld.tiff", height=7, width=9, units="in", res=300)
FeaturePlot(data_ld, features = c("SGK1"), label = FALSE) +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    plot.title = element_blank()  
  )
dev.off()

tiff("feature_PRKACA_data_ld.tiff", height=7, width=9, units="in", res=300)
FeaturePlot(data_ld, features = c("PRKACA"), label = FALSE) +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    plot.title = element_blank()  
  )
dev.off()

tiff("feature_PRKG1_data_ld.tiff", height=7, width=9, units="in", res=300)
FeaturePlot(data_ld, features = c("PRKG1"), label = FALSE) +
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    plot.title = element_blank()  
  )
dev.off()

# Calculate the percentage of cells expressing each gene among PT cells
# Extract only PT cells
pt_cells <- subset(data_ld, idents = "PT")

# Gene list
genes <- c("SGK1", "PRKACA", "PRKG1")

# Calculate the percentage of cells expressing each gene
expressing_percentages <- lapply(genes, function(gene) {
  # Get gene expression data
  expr_values <- FetchData(pt_cells, vars = gene)
  
  # Calculate the percentage of cells with expression greater than 0
  percent_expressing <- mean(expr_values[[gene]] > 0) * 100
  
  return(percent_expressing)
})

# Assign gene names as names to the list
names(expressing_percentages) <- genes

# Display results
print(expressing_percentages)