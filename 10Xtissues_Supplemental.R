#10X_Supplemental (tissues)

#Theme
library(randomcoloR)
#color_use <- distinctColorPalette(15) #length(levels(donor.integrated$subset_names))
color_use <- c("#7B52DD", "#7DE3D4", "#DC694C", "#77E595", "#E4DB61", "#D78AD6", "#CFDFDE", "#E0E2AC", "#9B9C59", "#778BD5", "#95EA4F", "#D3B6DB", "#75ABC5", "#DA4BA0", "#D49794")

#Directories
dir <- file.path("output")
out_dir <- file.path("figure images", "Supplemental")

#### Fig S2D (Integration donor)
DefaultAssay(donor.integrated) <- "integrated"

plot_umapdonor <- DimPlot(donor.integrated, group.by = "donor", pt.size = 0.5) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
file_name <- paste0("UMAP_origident.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_umapdonor, units = c("cm"))

plot_umapdonor <- DimPlot(donor.integrated, group.by = "donor", pt.size = 0.7, shape.by = "tissue") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none")
file_name <- paste0("UMAP_origident_tissuebyshape.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_umapdonor, units = c("cm"))


#### Tables (DE genes per subset)
idents_list <- levels(donor.integrated$subset_names)

de_pos_list <- {}
de_neg_list <- {}

for (i in 1:length(idents_list)){
  ident1 <- idents_list[i]
  
  condition.diffgenes <- FindMarkers(donor.integrated, ident.1 = ident1)
  #file_name <- paste0("cluster", ident1, "_degenes.csv")
  #write.csv(condition.diffgenes, file=file.path(out_dir, file_name))
  pos_dat <- subset(condition.diffgenes, avg_logFC > 0)
  neg_dat <- subset(condition.diffgenes, avg_logFC < 0)
  
  pos_marker <- rownames(pos_dat)
  neg_marker <- rownames(neg_dat)
  
  de_pos_list[[ident1]] <- pos_marker
  de_neg_list[[ident1]] <- neg_marker
}

de_pos_summary <- t(plyr::ldply(de_pos_list, rbind, .id = "cluster"))
de_neg_summary <- t(plyr::ldply(de_neg_list, rbind, .id = "cluster"))

file_name <- paste0("pos_degenes_summary.csv")
write.csv(de_pos_summary, file=file.path(out_dir, file_name), col.names = F)

file_name <- paste0("neg_degenes_summary.csv")
write.csv(de_neg_summary, file=file.path(out_dir, file_name), col.names = F)


for (i in 1:length(idents_list)){
  ident1 <- idents_list[i]
  
  condition.diffgenes <- FindMarkers(donor.integrated, ident.1 = ident1, only.pos = TRUE)
  
  file_name <- paste0("DEG_", ident1, ".csv")
  write.csv(condition.diffgenes, file=file.path(out_dir, file_name), col.names = F)
}

#### Fig S2F (DN differences)
DefaultAssay(donor.integrated) <- "RNA"

VlnPlot_2 <- function(object, features, idents, cols, pt.size, assay, xlab) {
  
  # Main function
  main_function <- function(object = object, features = features, idents = idents, cols = cols, pt.size = pt.size, assay = assay, xlab = xlab) {
    VlnPlot(object = object, features = features, idents = idents, cols = cols,  pt.size = pt.size, assay = assay) + 
      labs(x = xlab) + theme(legend.position = 'none')
  }
  
  # Apply main function on all features
  p <- lapply(X = features, object = object, idents = idents, cols = cols,  pt.size = pt.size, assay = assay, xlab = xlab,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument adapted from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  #cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(features))))
  cowplot::plot_grid(plotlist = p, ncol = 5)
}

dn_markers <- c("CXCR5", "ITGAX", "CR2", "CD19", "CD38", "CD24", "CD22", "CD69", "CD86", "HLA-DRA")

dn_topgenes <- VlnPlot_2(donor.integrated, 
                         features = dn_markers, 
                         idents = c("DN-A", "DN-B"), cols = c("#CFDFDE", "#D78AD6"), pt.size = 0.2, assay = "RNA", xlab = " ")

file_name <- paste0("dn_vlnplot_select.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = dn_topgenes, width = 12, height = 8)

## DotPlot for aNAV
DefaultAssay(donor.integrated) <- "RNA"

anav_markers <- c("CD19", "CXCR5", "CD24", "CD38", "ITGAX")
anav_idents <- c("aNAV", "Naive", "TS")

plot_anavdot <- DotPlot(donor.integrated, features = anav_markers, scale.by = "radius", idents = anav_idents, dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  scale_size(range = c(2, 10)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_anav_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_anavdot, width = 12, height = 12, units = c("cm"))
