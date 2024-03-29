#Figure2
#Code for transcriptomic analysis shown in Figure 2 (a-d, f) and Fig 6 (h, j)

#Library
library(tidyr)
library(Seurat) #V3.2.2
library(patchwork)
library(tidyverse)

library(textTinyR) #calculate sparse matrix stats
library(rstatix)
library(dunn.test)
library(car)

library(ggplot2)
library(viridis)
library(ggpubr)

#Theme
library(randomcoloR)
color_use <- c("#7B52DD", "#7DE3D4", "#DC694C", "#77E595", "#E4DB61", "#D78AD6", "#CFDFDE", "#E0E2AC", "#9B9C59", "#778BD5", "#95EA4F", "#D3B6DB", "#75ABC5", "#DA4BA0", "#D49794")

#Directories
dir <- file.path("R_analysis", "output")
dir_clust <- file.path("Data", "Seurat saved objects")

out_dir <- file.path("figure images", "Fig2")

#Import data
donor.integrated <- readRDS(file = file.path(dir, "integrated_tissue.rds"))

#Clustering
pca_dim <- 1:25
cluster_res <- 1.3

DefaultAssay(donor.integrated) <- "integrated"

donor.integrated <- RunPCA(donor.integrated)
donor.integrated <- RunUMAP(object = donor.integrated, dims = pca_dim, n.neighbors = 30)
donor.integrated <- FindNeighbors(donor.integrated, reduction = "pca", dims = pca_dim, 
                                  assay = "integrated", nn.method="rann") #nn.method="rann"
donor.integrated <- FindClusters(donor.integrated, resolution = cluster_res)

##Rearrange order based on cluster tree (hierarchical)
donor.integrated <- BuildClusterTree(object = donor.integrated, dims = pca_dim)

data.tree <- Tool(donor.integrated, slot = "BuildClusterTree")
is_tip <- data.tree$edge[,2] <= length(data.tree$tip.label)
ordered_tips <- data.tree$edge[is_tip, 2]

my_levels <- data.tree$tip.label[ordered_tips]
donor.integrated@active.ident <- factor(x = donor.integrated@active.ident, levels = my_levels)

##RNA scaling
donor.integrated <- NormalizeData(donor.integrated, assay = "RNA")
all.genes <- rownames(donor.integrated)
donor.integrated <- ScaleData(donor.integrated, features = all.genes, split.by = "donor", assay = "RNA")

##ADT Scaling
donor.integrated <- NormalizeData(donor.integrated, assay = "ADT", normalization.method = "CLR")
donor.integrated <- ScaleData(donor.integrated, assay = "ADT", split.by = "orig.ident")

DefaultAssay(donor.integrated) <- "RNA"
donor.integrated <- subset(donor.integrated, subset = CD14 < 1.2)

#Rename subsets
new.cluster.ids <- c("DN-A", "ABC3", "DN-B", "ABC1", "GC", "GC", "GC", "IgM-only", "Memory", "Memory", "MZB-1", "Memory", "TS", "Naive", "Naive", "aNAV", "IgM-only", "Naive", "ABC2", "aNAV", "MZB-2", "ABC4", "PB")
names(new.cluster.ids) <- levels(donor.integrated)
donor.integrated <- RenameIdents(donor.integrated, new.cluster.ids)

donor.integrated[["subset_names"]] <- donor.integrated@active.ident
donor.integrated@active.ident <- factor(donor.integrated@active.ident, levels = sort(levels(donor.integrated$subset_names), decreasing = T)) #add , decreasing = T to sort
donor.integrated$subset_names <- factor(donor.integrated$subset_names, levels = sort(levels(donor.integrated$subset_names)))

saveRDS(donor.integrated, file = file.path(out_dir, "donorintegrated.RDS"))

####Fig2a (ADT vs RNA UMAP)
adt_markers <- c("adt_CD27", "adt_CD38", "adt_IGM", "adt_IGD")
rna_markers <- c("CD27", "CD38", "IGHM", "IGHD")

#ADT
cite_plot <- FeaturePlot(donor.integrated, features = adt_markers, min.cutoff = "q10", max.cutoff = "q95", ncol = 1, cols = viridis(200), order = T, combine = TRUE) & NoAxes() & NoLegend()

file_name <- paste0("umap_citeseqmarkers.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = cite_plot, width = 10, height = 40, units = c("cm"))

#RNA
rna_plot <- FeaturePlot(donor.integrated, features = rna_markers, min.cutoff = "q15", max.cutoff = "q95", ncol = 1, cols = viridis(200), order = T) & NoAxes() & NoLegend()

file_name <- paste0("umap_rnaseqmarkers.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = rna_plot, width = 10, height = 40, units = c("cm"))


####Fig2b (subset key)
DefaultAssay(donor.integrated) <- "integrated"

plot_umaptissue <- DimPlot(donor.integrated, reduction = "umap", label = TRUE, pt.size = 0.7, cols = color_use) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
file_name <- paste0("UMAP_subsetkey.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_umaptissue, units = c("cm"))

####Fig2c (dot plot)
DefaultAssay(donor.integrated) <- "RNA"

donor.integrated <- AddModuleScore(donor.integrated, list(c("IGHA1", "IGHA2")), name='IgA.RNA', assay = "RNA")
donor.integrated <- AddModuleScore(donor.integrated, list(c("IGHG1", "IGHG2", "IGHG3", "IGHG4")), name='IgG.RNA', assay = "RNA")
donor.integrated <- AddModuleScore(donor.integrated, list(c('MME', 'BCL6', 'BCL7A', "PCNA")), name='GC.RNA', assay = "RNA")

rna_features_output_split_dot <- c("IGHD", "IGHM", "CD27", "CD1C", "CD24", "CD38", "CR2", "BCL6", "BCL7A", "PCNA", "MME", "MZB1", "COCH", "PLD4", "ZEB2", "CCR7", "CXCR5", "IGLL5", "MX1", "GC.RNA1", "IgG.RNA1", "IgA.RNA1")
cytof_markers <- c("CD180", "MS4A1", "CD40", "IGHM", "PTPRC", "CD27", "CD24", "TNFSF13B", "CXCR5", "CCR7", "CD38", "IL2RA", "IgA.RNA1", "IGHD", "TNFRSF17", "TNFRSF13B", "MME", "IL7R", "CD80", "PDCD1", "ITGB7")
cytof_markers2 <- c("CCR7", "CXCR5", "TNFRSF13C", "CD24", "CD27")
other <- c("DDX21")

plot_cytofdot <- DotPlot(donor.integrated, features = cytof_markers, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))
  
file_name <- paste0("tissues_cytofdotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 21, height = 15, units = c("cm"))

plot_cytofdot <- DotPlot(donor.integrated, features = cytof_markers2, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_sig_cytofdotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 12, height = 16, units = c("cm"))

#Fig 6J
#b7 only
plot_cytofdot <- DotPlot(donor.integrated, features = "ITGB7", scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_b7_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 7, height = 16, units = c("cm"))

plot_cytofdot <- DotPlot(donor.integrated, features = "B7", scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", scale.min = 80, scale.max = 100) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_discrete(labels = "\u03B27_adt") +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_adtb7_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 10, height = 20, units = c("cm"))

plot_cytofdot <- DotPlot(donor.integrated, features = c("ITGB7", "B7"), scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", scale.min = 25, scale.max = 100) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_discrete(labels = c("ITGB7", "\u03B27_adt")) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_b7_rnaadt_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 15, height = 20, units = c("cm"))

#
rna_features_output_split_dot <- c("IGHD", "IGHM", "CD27", "CD1C", "PLD4", "CR2", "CCR7", "MME", "MZB1", "IGLL5", "CD24", "CD38", "COCH", "CXCR5", "ZEB2", "MX1", "BCL6", "BCL7A", "PCNA", "GC.RNA1", "IgG.RNA1", "IgA.RNA1")


plot_alldot <- DotPlot(donor.integrated, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_alldot, width = 25, height = 15, units = c("cm"))

plot_alldot <- DotPlot(donor.integrated, split.by = "tissue", cols = c("blue", "black", "red"),
                       features = other, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
file_name <- paste0("tissues_dotplot_DDX21.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_alldot, width = 25, height = 15, units = c("cm"))


####Fig2d (UMAP by tissue)
DefaultAssay(donor.integrated) <- "integrated"

plot_umaptissue <- DimPlot(donor.integrated, reduction = "umap", label = TRUE, pt.size = 0.7, cols = color_use, split.by = "tissue") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
file_name <- paste0("UMAP_subset_bytissue.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_umaptissue, width = 35, units = c("cm"))

####Fig2f (Abundance)
meta.data <- donor.integrated@meta.data

counts <- meta.data %>% group_by(tissue, donor, subset_names) %>% summarise(count = n())
percentage <- counts %>% group_by(tissue, donor) %>% mutate(percent = count/sum(count)*100)

bxp <- ggboxplot(percentage, x = "subset_names", y = "percent",
  color = "tissue", palette = c("APP" = "blue","MLN" = "black","SPL" = "red"), add = "jitter", shape = "tissue", outlier.shape = NA
) + labs(x = "Subset Names", y = "% of CD19", color = "Tissue", shape = "Tissue") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title.x=element_blank())

file_name <- paste0("Abundance_stats.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = bxp, width = 25, height = 9, units = c("cm"))

#Stats
cd19.model <- lm(percent ~ subset_names*tissue, data = percentage)
percentage %>%
  group_by(subset_names) %>%
  anova_test(percent ~ tissue, error = cd19.model)

pwc <- percentage %>% 
  group_by(subset_names) %>%
  emmeans_test(percent ~ tissue, p.adjust.method = "bonferroni") 
pwc

#Fig 6H
donor.integrated@active.ident <- factor(donor.integrated@active.ident, levels = sort(levels(donor.integrated$subset_names), decreasing = T)) #add , decreasing = T to sort

cytof_markers <- c("CCR7", "CXCR5", "TNFRSF13C", "CD24", "CD27")

plot_cytofdot <- DotPlot(donor.integrated, features = cytof_markers, scale.by = "radius", dot.min = 0, dot.scale = 15, col.min = -2, col.max = 2, assay = "RNA", scale.min = 0, scale.max = 60) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_sig_cytofdotplot_allsame.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 12, height = 16, units = c("cm"))


#Fig 6J
plot_cytofdot <- DotPlot(donor.integrated, features = "B7", scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", scale.min = 80, scale.max = 100) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_x_discrete(labels = "\u03B27_adt") +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))

file_name <- paste0("tissues_adtb7_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 10, height = 20, units = c("cm"))
