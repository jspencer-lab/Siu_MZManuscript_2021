#Figure 4

#Library
library(tidyr)
library(Seurat)
library(patchwork)
library(tidyverse)

library(ggplot2)
library(viridis)
library(ggpubr)

#Directories
dir <- file.path("")
out_dir <- file.path("figure images", "Fig4")

#Import data
donor.integrated <- readRDS(file = file.path(dir, "integration_allsamples", "integrated_tissue.rds"))

#### Fig 4a (volcano)
library(EnhancedVolcano)
mz.de.markers <- FindMarkers(donor.integrated, ident.1 = "MZB-2", ident.2 = "MZB-1", logfc.threshold = 0)

markerslabel <- c("DDX21", "TXNIP", "MIF", "NCL", "NME1", "FOS", "HLA-DRB1", "YBX1", "HLA-DQA1", "PSME2", "HLA-DRA", "PTMA", "FABP5", "HSP90AB1",
                  "TXNIP", "MALAT1", "LINC00926", "LTB", "CD37", "CD27", "AIM2", "CD83")

mz_volcano <- EnhancedVolcano(mz.de.markers, lab = rownames(mz.de.markers), 
                x = "avg_logFC", y = "p_val",
                pCutoff = 10e-6, FCcutoff = 0.25, 
                title = 'MZB-1 versus MZB-2', subtitle = NULL,
                labSize = 6, 
                selectLab = markerslabel, 
                labCol = 'black', labFace = 'bold',
                boxedLabels = T, drawConnectors = T)
file_name <- paste0("mzb1vs2_volcano_select.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = mz_volcano, width = 10, height = 10)

#### Fig 4b (vln select)
VlnPlot_2 <- function(object, features, idents, split.by, cols, pt.size, assay, xlab) {
  
  # Main function
  main_function <- function(object = object, features = features, idents = idents, split.by = split.by, cols = cols, pt.size = pt.size, assay = assay, xlab = xlab) {
    VlnPlot(object = object, features = features, idents = idents, split.by = split.by, cols = cols,  pt.size = pt.size, assay = assay) + 
      labs(x = xlab) + theme(legend.position = 'none')
  }
  
  # Apply main function on all features
  p <- lapply(X = features, object = object, idents = idents, split.by = split.by, cols = cols,  pt.size = pt.size, assay = assay, xlab = xlab,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument adapted from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(features))))
}

mz_topgenes <- VlnPlot_2(donor.integrated, 
                       features = c("DDX21", "CD83", "MIF", "TXNIP", "LTB", "CD37"), 
                       idents = c("MZB-1", "MZB-2"), split.by = "tissue", cols = c("blue", "black", "red"), pt.size = 0.2, assay = "RNA", xlab = " ")

file_name <- paste0("mz_vlnplot_select.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = mz_topgenes, width = 12, height = 8)

#Fig 4C was created using GSEA software from Broad Institute (https://www.gsea-msigdb.org/gsea/index.jsp)

