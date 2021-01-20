library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(RcmdrMisc) #rcorr.adjust function
library(corrplot) #to make correlation plots
library(ggpubr)


# Load -------------------------------------------------------------
# Data
integrated_samples  <- readRDS("blood_integrated_data.rds")

# Global vars
SAMPLES <- c("health_1", "health_2", "health_3", "lupus_1", 'lupus_2', 'lupus_3')

MANUAL_CLASSIFICATION = 'manual_classification'
CLASSIFICATIONS <- list()
CLASSIFICATIONS[['TRANSITIONAL']] = 'TS'
CLASSIFICATIONS[['NAIVE']] = 'Naive'
CLASSIFICATIONS[['MZB_2']] = 'MZB-2'
CLASSIFICATIONS[['MZB_1']] = 'MZB-1'
CLASSIFICATIONS[['IGM_ONLY']] = 'IgM-only'
CLASSIFICATIONS[['CSM']] = 'Memory'
CLASSIFICATIONS[['ACTIVATED_NAIVE']] = 'aNAV'
CLASSIFICATIONS[['DN2']] = 'DN'
CLASSIFICATIONS[['PLASMABLASTS']] = 'PB'
CLASSIFICATIONS[['MALAT1']] = 'ABC3'

# Format
COLOURS <- list()
COLOURS[[MANUAL_CLASSIFICATION]] <- list()
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$TRANSITIONAL]] =  "#D49794"
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$NAIVE]] = "#75ABC5"
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$MZB_2]] = "#D3B6DB"
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$MZB_1]] = "#95EA4F"
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$IGM_ONLY]] ="#9B9C59" 
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$CSM]] = '#778BD5'
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$ACTIVATED_NAIVE]] = "#E4DB61" 
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$DN2]] = "#CFDFDE"
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$PLASMABLASTS]] = "#DA4BA0" 
COLOURS[[MANUAL_CLASSIFICATION]][[CLASSIFICATIONS$MALAT1]] = "#DC694C"
COLOURS[['status']] <- c("health" = "darkgreen", "lupus" = "orange")

# Dir
out_dir <- 'plots/Figure6'

# UMAP with classification -------------------------------------------------------
panel_6a <- DimPlot(integrated_samples, group.by = MANUAL_CLASSIFICATION, cols = as.character(COLOURS$manual_classification[sort(names(COLOURS$manual_classification))]), label = TRUE, label.size = 8, pt.size=0.7) + theme(axis.ticks = element_blank(), axis.text = element_blank())
file_name <- ("Blood_UMAP_with_classification.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = panel_6a, width = 35, units = c("cm"))

# CD1C and PLD4 -------------------------------------------------------
panel_6b <- (FeaturePlot(integrated_samples, features = c("PLD4"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, cols = viridis(200), order = T, pt.size = 1.7) & NoAxes() & NoLegend()) | (FeaturePlot(integrated_samples, features = c("CD1C"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, cols = viridis(200), order = T, pt.size = 1.7) & NoAxes() + theme(legend.text = element_blank()))
file_name <- ("PLD4_CD1C.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = panel_6b, width = 60, units = c("cm"))

# Mu freq -------------------------------------------------------
plot_data <- integrated_samples@meta.data[!is.na(integrated_samples$mu_freq), c('mu_freq', MANUAL_CLASSIFICATION)]
panel_6d <- ggplot(plot_data, aes(x=manual_classification, y=mu_freq, fill=manual_classification)) + 
            geom_boxplot() + 
            scale_x_discrete(limits=sort(as.character(CLASSIFICATIONS))) + 
            ylab("Mutation frequency") + 
            scale_fill_manual(values=COLOURS$manual_classification) + 
            theme_classic() + 
            theme(axis.title.x = element_blank(), text = element_text(size=20), axis.text.x=element_text(angle=45,hjust=1), 
              legend.position = "none")
file_name <- ("mutation_frequency_by_classification.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = panel_6d, width = 35, units = c("cm"))

# clonal correlations -------------------------------------------
metadata <- integrated_samples@meta.data
metadata$manual_classification <- as.factor(metadata$manual_classification)

classif_order <- sort(as.character(CLASSIFICATIONS))

clone_tables <- list()
for (sample in SAMPLES){
  g <- metadata[metadata$Sample == sample & !(is.na(metadata$clone_id_HC)),]
  g <- g[g$clone_id_HC %in% names(which(table(g$clone_id_HC) > 1)),]
  g <- g[,c('clone_id_HC', MANUAL_CLASSIFICATION)] 
  # cross_classification_clones <- names(which(rowSums(table(g) > 0) > 1))
  # g <- g[g$clone_id %in% cross_classification_clones,]
  clone_tables[[sample]] <- as.data.frame.matrix(table(g))
  rownames(clone_tables[[sample]]) <- paste0(sample, '_', rownames(clone_tables[[sample]]))
  for (c in colnames(clone_tables[[sample]])){
    clone_tables[[sample]][[c]] <- as.integer(clone_tables[[sample]][[c]] > 0)
  }
}

# HEALTHY CORRELATION
healthy <- rbind(clone_tables$health_1, clone_tables$health_2, clone_tables$health_3)
healthy_cor <- rcorr.adjust(as.matrix(healthy), type="spearman")

healthy_cor$P[healthy_cor$P == "<.0001"] <- 0.0001
healthy_cor$P<- data.matrix(healthy_cor$P)
class(healthy_cor$P) <- "numeric"

healthy_cor$R$r <- healthy_cor$R$r[classif_order, classif_order]
healthy_cor$R$P <- healthy_cor$R$P[classif_order, classif_order]
healthy_cor$P  <- healthy_cor$P[classif_order, classif_order]

pdf(file.path(out_dir, "correlation_health.pdf"), width = 10)
corrplot(healthy_cor$R$r, type="lower", order="original", method = "color", 
         p.mat = healthy_cor$P, sig.level = c(.001, .01, .05), pch.cex = 1, 
         insig = "label_sig", tl.col = "black", tl.srt = 90, outline = "gray",
         title = paste0("HEALTHY"), mar=c(0,0,1,0), diag = FALSE)
dev.off()


# LUPUS CORRELATION
lupus <- rbind(clone_tables$lupus_1, clone_tables$lupus_2, clone_tables$lupus_3)
lupus_cor <- rcorr.adjust(as.matrix(lupus), type="spearman")

lupus_cor$P[lupus_cor$P == "<.0001"] <- 0.0001
lupus_cor$P<- data.matrix(lupus_cor$P)
class(lupus_cor$P) <- "numeric"

lupus_cor$R$r <- lupus_cor$R$r[classif_order, classif_order]
lupus_cor$R$P <- lupus_cor$R$P[classif_order, classif_order]
lupus_cor$P  <-  lupus_cor$P[classif_order, classif_order]

pdf(file.path(out_dir, "correlation_lupus.pdf"), width = 10)
corrplot(lupus_cor$R$r, type="lower", order="original", method = "color", 
         p.mat = lupus_cor$P, sig.level = c(.001, .01, .05), pch.cex = 1, 
         insig = "label_sig", tl.col = "black", tl.srt = 90, outline = "gray",
         title = paste0("HEALTHY"), mar=c(0,0,1,0), diag = FALSE)
dev.off()

# Dotplots 2 ----
temp <- subset(integrated_samples, cells = WhichCells(integrated_samples, expression = status == 'health'))

temp$manual_classification <- factor(temp$manual_classification, levels = sort(as.character(CLASSIFICATIONS), decreasing = T))
plot_cytofdot <- DotPlot(temp, features = "ITGB7", scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", 
                         group.by = MANUAL_CLASSIFICATION) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))
file_name <- paste0("blood_b7_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 7, height = 16, units = c("cm"))

cytof_markers2 <- c("CCR7", "CXCR5", "TNFRSF13C", "CD24", "CD27")
plot_cytofdot <- DotPlot(temp, features = cytof_markers2, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", group.by = MANUAL_CLASSIFICATION) +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"))
file_name <- paste0("blood_sig_cytofdotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = plot_cytofdot, width = 12, height = 16, units = c("cm"))

DefaultAssay(temp) <- 'ADT'
adtb7 <- DotPlot(temp, features = c("ADT-B7-TotalSeqC"), scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", 
                 group.by = MANUAL_CLASSIFICATION) +
  scale_color_viridis_c() +
  RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression")) +  
  scale_x_discrete(label=c(expression(beta*"7_adt")))
file_name <- paste0("blood_b7_adt_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = adtb7, width = 10, height = 16, units = c("cm"))

DefaultAssay(temp) <- 'RNA'
b7 <- DotPlot(temp, features = c("ITGB7", "ADT-B7-TotalSeqC"), scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA", 
                         group.by = MANUAL_CLASSIFICATION) +
  scale_color_viridis_c() +
  RotatedAxis() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression")) +
  scale_x_discrete(label=c('ITGB7', expression(beta*"7_adt")))
file_name <- paste0("blood_b7_rnaadt_dotplot.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = b7, width = 10, height = 16, units = c("cm"))

# IgM ---------------------------------------------

panel_6g <- FeaturePlot(integrated_samples, features = c("ADT-IGM-TotalSeqC"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, cols = viridis(200), pt.size = 0.7, order=T, cells = WhichCells(integrated_samples, expression = status == 'health')) & NoLegend() & ggtitle("adt_IGM") & theme(axis.ticks = element_blank(), axis.text = element_blank())
file_name <- ("adt_IGM_health.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = panel_6g, width = 30, units = c("cm"))

# Abundance ---------------------------------------------
meta.data <- integrated_samples@meta.data

counts <- meta.data %>% group_by(status, Sample, manual_classification) %>% summarise(count = n())
percentage <- counts %>% group_by(status, Sample) %>% mutate(percent = count/sum(count)*100)

pdf(file.path(out_dir, "abundances.pdf"), width = 7, height = 5)
ggboxplot(percentage, x = MANUAL_CLASSIFICATION, y = "percent",
                 color = "status", palette = COLOURS$status, add = "jitter", shape = "status",
                 outlier.shape = NA) + labs(x = "Subset Names", y = "% of CD19", color = "Status", shape = "Status") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        axis.title.x=element_blank())
dev.off()

#Stats
cd19.model <- lm(percent ~ manual_classification*status, data = percentage)
percentage %>%
  group_by(manual_classification) %>%
  anova_test(percent ~ status, error = cd19.model)

pwc <- percentage %>% 
  group_by(manual_classification) %>%
  emmeans_test(percent ~ status, p.adjust.method = "bonferroni") 
pwc

x = percentage$percent[percentage$status == 'health' & percentage$manual_classification == 'MZB-1']
y = percentage$percent[percentage$status == 'lupus' & percentage$manual_classification == 'MZB-1']

wilcox.test(x = x, y=y)
