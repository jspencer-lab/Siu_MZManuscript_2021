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
dir <- file.path("C:", "Users", "jacqu", "OneDrive - University Of Cambridge", "PhD", "10X", "R_analysis", "output")
out_dir <- file.path("C:", "Users", "jacqu", "Dropbox", "NI Submission 2021", "figure images", "Fig4")

#Import data
donor.integrated <- readRDS(file = file.path(dir, "JS10X_012 (Integration with regression)", "integration_allsamples", "integrated_tissue.rds"))

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

#### Fig 4a (linear)
mz.cells <- subset(donor.integrated, idents = c("MZB-1", "MZB-2"))

avg.mz.cells <- log1p(AverageExpression(mz.cells, verbose = FALSE)$RNA)
avg.mz.cells$gene <- rownames(avg.mz.cells)
avg.mz.cells$diff <- abs(avg.mz.cells$`MZB-1` - avg.mz.cells$`MZB-2`)
avg.mz.cells$`MZB-1` <- as.numeric(avg.mz.cells$`MZB-1`)
avg.mz.cells$`MZB-2` <- as.numeric(avg.mz.cells$`MZB-2`)

genesdiff <- subset(avg.mz.cells, diff > 0.5 & gene == "CCR7")

genesdiff <- subset(avg.mz.cells, gene == "CCR7")


p1 <- ggplot(avg.mz.cells, aes(x = `MZB-1`, y = `MZB-2`)) + 
      geom_point() + 
      geom_smooth(method=lm , color="navy", fill="grey", se=TRUE) +
      theme_pubr()
p1 <- LabelPoints(plot = p1, points = genesdiff$gene, repel = TRUE, xnudge = 0, ynudge = 0, colour = "red")

file_name <- paste0("MZB1vs2_linear.png")
ggsave(filename = file.path(out_dir, file_name), plot = p1, width = 30, height = 30, units = c("cm"))

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

#### Fig 4x (mz dotplot)
dot_marker <- c("DDX21", "CD83", "MIF", "NCL", "TNFRSF1B", "CCR7", "CD24", "CR2", "TXNIP", "LTB", "CD37", "CD27")

plot_mzdot <- DotPlot(donor.integrated, features = dot_marker, idents = c("MZB-2", "MZB-1"), scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_size(range = c(1, 10)) +
  #scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
file_name <- paste0("tissues_mzdotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_mzdot, width = 20, height = 12, units = c("cm"))

#### Fig4x (tbet dotplot)
tbet_dot <- DotPlot(donor.integrated, features = "TBX21", scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", group.by = "subset_names", split.by = "tissue", cols = c("blue", "black", "red"), 
                    idents = c("MZB-1", "MZB-2", "aNAV", "ABC1", "ABC2", "DN-A", "DN-B", "Memory", "GC")) + labs(x = " ", y = " ") 
file_name <- paste0("tbet_dotplot_short.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = tbet_dot, width = 5, height = 10)

#V2
tbet_dot_dat <- tbet_dot$data
id <- tbet_dot_dat$id
tbet_dot_dat$subset <- unlist(stringr::str_split(id, "_"))[ c(TRUE,FALSE) ]
tbet_dot_dat$tissue <- unlist(stringr::str_split(id, "_"))[ c(FALSE,TRUE) ]

mid <- mean(tbet_dot_dat$avg.exp)
avg.exp.lim <- tbet_dot_dat$avg.exp
avg.exp.lim[avg.exp.lim > 0.1] <- 0.1
tbet_dot_dat$avg.exp.lim <- avg.exp.lim
    
tbet_dot_V2 <- tbet_dot_dat %>% 
  ggplot(aes(x=tissue, y = subset, color = avg.exp.lim, size = pct.exp*2)) + 
  scale_color_viridis(option = "D", limits = c(0,0.1)) +
  geom_point() +
  theme_pubr() +
  labs(colour = "Average \nTbet \nExpression", size = "% Expressed", x = " ", y = " ") +
  theme(legend.position = "right")
file_name <- paste0("tbet_dotplot_maxed.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = tbet_dot_V2, width = 10, height = 10, units = c("cm"))

tbet_dot_V2 <- tbet_dot_dat %>% 
  ggplot(aes(x=tissue, y = subset, color = avg.exp, size = pct.exp*2)) + 
  scale_color_viridis(option = "D") +
  geom_point() +
  theme_pubr()
file_name <- paste0("tbet_dotplot_V2.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = tbet_dot_V2, width = 10, height = 10, units = c("cm"))

#### Fig4x (fcrl4 dotplot)
fcrl4_dot <- DotPlot(donor.integrated, features = "FCRL4", scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", group.by = "subset_names", split.by = "tissue", cols = c("blue", "black", "red"), 
                    idents = c("MZB-1", "MZB-2", "aNAV", "ABC1", "ABC2", "DN-A", "DN-B", "Memory", "GC")
                    ) + labs(x = " ", y = " ") 
file_name <- paste0("fcrl4_dotplot_short.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = fcrl4_dot, width = 5, height = 10)

#V2
dot_dat <- fcrl4_dot$data
id <- dot_dat$id
dot_dat$subset <- unlist(stringr::str_split(id, "_"))[ c(TRUE,FALSE) ]
dot_dat$tissue <- unlist(stringr::str_split(id, "_"))[ c(FALSE,TRUE) ]

fcrl4_dot_V2 <- dot_dat %>% 
  ggplot(aes(x=tissue, y = subset, color = avg.exp, size = pct.exp*2)) + 
  scale_color_viridis(option = "D") +
  geom_point() +
  theme_pubr() +
  labs(colour = "Average \nFcRL4 \nExpression", size = "% Expressed", x = " ", y = " ") +
  theme(legend.position = "right")
file_name <- paste0("fcrl4_dotplot_V2.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = fcrl4_dot_V2, width = 10, height = 11, units = c("cm"))