library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)


# Load -------------------------------------------------------------
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
COLOURS[['status']] <- c("health" = "dodgerblue", "lupus" = "orange")

# Dir
out_dir <- 'plots/Supp'

# Load 10X Data  -------------------------------
# Load the 10X data from each folder and create a Seurat object

health_samples = list()
lupus_samples = list()

# Minimum number of cells a gene must be present in to be included
QC_MIN_GENES_PER_CELL <- 3
# Genes denoting B Cells
BCELL_GENES <- c("CD79A", "CD79B", "CD19", "MS4A1")
# Variable region genes (which will be excluded - moved to a separate assay)
IG_GENES <- c("^IGHV", "^IGLV", "^IGLC", "^IGKV", "^IGKC", "^AC233755")

for (sample in SAMPLES){
  # Read the 10X data, should be 2 matrices - Gene Expression and Antibody Capture
  temp_data <- Read10X(paste('data/raw/', sample, '/outs/filtered_feature_bc_matrix/', sep=''))
  
  # Find the rows which are IG genes
  ig_rows = c()
  for (pattern in IG_GENES){
    ig_rows <- c(ig_rows, grep(pattern, x = rownames(temp_data$`Gene Expression`)))
  }
  
  # Remove the IG gene related rows and split them out into a new assay
  temp_data$`IG genes` <- temp_data$`Gene Expression`[ig_rows,]
  temp_data$`Gene Expression` <- temp_data$`Gene Expression`[-ig_rows,]
  
  # Create the seurat object
  temp_seurat <- CreateSeuratObject(counts = temp_data[['Gene Expression']], project = sample, min.cells = QC_MIN_GENES_PER_CELL)
  
  # Add antibody capture
  temp_seurat[["ADT"]] <- CreateAssayObject(temp_data[["Antibody Capture"]][, colnames(x = temp_seurat)])
  
  # Add IG genes - this will sit in a separate assay, and won't impact UMAP etc. but can be used for FeaturePlots
  temp_seurat[["IG"]] <- CreateAssayObject(temp_data[["IG genes"]][, colnames(x = temp_seurat)])
  
  # Add the sample name as a meta data column (for identification downstream when data is integrated)
  temp_seurat[['sample']] <- sample
  
  # Save the barcode in the meta.data and add the sample ID to the cell name (for downstream integration)
  temp_seurat[['barcode']] <- colnames(temp_seurat)
  temp_seurat <- RenameCells(object = temp_seurat, add.cell.id = sample)
  
  # Add mitochondrial percent
  temp_seurat[["percent.mt"]] <- PercentageFeatureSet(temp_seurat, pattern = "^MT-")
  # Add ribosomal percent
  temp_seurat[["percent.rb"]] <- PercentageFeatureSet(temp_seurat, "^RP[SL]")
  # Add B cell percent
  temp_seurat[["percent.Bcell"]] <- PercentageFeatureSet(object = temp_seurat, features = BCELL_GENES)
  
  # Add to appropriate list and add a status value based on sample name
  if (substr(sample, 1, 1) == 'h'){
    temp_seurat[['status']] <- 'health'
    health_samples[[sample]] <- temp_seurat
  } else {
    temp_seurat[['status']] <- 'lupus'
    lupus_samples[[sample]] <- temp_seurat
  }
}

# Classify CD19 / TS  -----------------------------------------------------
# Data is a mix of full CD19 and sorted transitional cells. Cells with AHH1 hash
# should be full CD19, cells with AHH2 hash should be sorted TS. So split the
# dataset based on expression of these.

# Thresholds
AHH1_thresholds_health <- list('health_1'= 30, 'health_2' = 100, 'health_3' = 120)
AHH2_thresholds_health <- list('health_1'= 3, 'health_2' = 30, 'health_3' = 25)
AHH1_thresholds_lupus <- list('lupus_1'= 30, 'lupus_2' = 150, 'lupus_3' = 120)
AHH2_thresholds_lupus <- list('lupus_1'= 3, 'lupus_2' = 30, 'lupus_3' = 25)

# Subset
health_CD19_samples <- list()
health_TS_samples <- list()
lupus_CD19_samples <- list()
lupus_TS_samples <- list()

# Split for health samples
for (sample in ls(health_samples)){
  temp_data <- health_samples[[sample]]
  DefaultAssay(temp_data) <- 'ADT'
  CD19_cell_names <- WhichCells(temp_data, expression = `HTO-AHH1-TotalSeqC` > AHH1_thresholds_health[[sample]] &
                                  `HTO-AHH2-TotalSeqC` <= AHH2_thresholds_health[[sample]])
  TS_cell_names <- WhichCells(temp_data, expression = `HTO-AHH1-TotalSeqC` <= AHH1_thresholds_health[[sample]] &
                                `HTO-AHH2-TotalSeqC` > AHH2_thresholds_health[[sample]])
  health_CD19_samples[[sample]] <- subset(temp_data, cells=CD19_cell_names)
  health_TS_samples[[sample]] <- subset(temp_data, cells=TS_cell_names)
}

# Repeat for lupus samples
for (sample in ls(lupus_samples)){
  temp_data <- lupus_samples[[sample]]
  DefaultAssay(temp_data) <- 'ADT'
  CD19_cell_names <- WhichCells(temp_data, expression = `HTO-AHH1-TotalSeqC` > AHH1_thresholds_lupus[[sample]] &
                                  `HTO-AHH2-TotalSeqC` <= AHH2_thresholds_lupus[[sample]])
  TS_cell_names <- WhichCells(temp_data, expression = `HTO-AHH1-TotalSeqC` <= AHH1_thresholds_lupus[[sample]] &
                                `HTO-AHH2-TotalSeqC` > AHH2_thresholds_lupus[[sample]])
  lupus_CD19_samples[[sample]] <- subset(temp_data, cells=CD19_cell_names)
  lupus_TS_samples[[sample]] <- subset(temp_data, cells=TS_cell_names)
}

# QUALITY CONTROL - FEATURES ----------------------------------------------
# Remove features that are no longer required (hashes) and store as metadata.
health_CD19_QC_samples <- list()
lupus_CD19_QC_samples <- list()

QC_features <- function(seurat_obj){
  # Get the hash data and save it into the meta.data
  hash_data <- FetchData(seurat_obj, c('HTO-AHH1-TotalSeqC', 'HTO-AHH2-TotalSeqC'))
  seurat_obj[['HASH1']] <- hash_data$`HTO-AHH1-TotalSeqC`
  seurat_obj[['HASH2']] <- hash_data$`HTO-AHH2-TotalSeqC`
  
  # Get a list of all features and remove the hashes from it
  features_to_keep <- c(rownames(seurat_obj@assays$ADT),
                        rownames(seurat_obj@assays$RNA),
                        rownames(seurat_obj@assays$IG))
  features_to_keep <- features_to_keep[features_to_keep != 'HTO-AHH1-TotalSeqC' & features_to_keep != 'HTO-AHH2-TotalSeqC']
  
  # Subset to just the kept features (i.e. remove the hashes from the ADT assay)
  return(subset(seurat_obj, features = features_to_keep))
}

# Health
for (sample in ls(health_CD19_samples)){
  health_CD19_QC_samples[[sample]] <- QC_features(health_CD19_samples[[sample]])
}

# Lupus
for (sample in ls(lupus_CD19_samples)){
  lupus_CD19_QC_samples[[sample]] <- QC_features(lupus_CD19_samples[[sample]])
}


# QUALITY CONTROL - CELLS -------------------------------------------------
# Perform quality control on cells. We QC for:
# Number of reads
# Number of genes
# % of reads which are mitochondrial
# % of reads which are B cell genes

# Function to automatically calculate QC thresholds by MAD (see ?mad)
QC_sample <- function(seurat_obj, sample_name, colour, filename){
  # nFeature_RNA
  v <- seurat_obj[['nFeature_RNA']]$nFeature_RNA
  v <- log(v)
  nFeature_min <- exp(median(v) - (3*mad(v)))
  nFeature_max <- exp(median(v) + (3*mad(v)))
  
  # nCount_RNA
  v <- seurat_obj[['nCount_RNA']]$nCount_RNA
  v <- log(v)
  nCount_min <- exp(median(v) - (3*mad(v)))
  nCount_max <- exp(median(v) + (3*mad(v)))
  
  # percent.mt
  v <- seurat_obj[['percent.mt']]$percent.mt
  v <- log(v)
  mt_max <- exp(median(v) + (3*mad(v)))
  
  # Bcell min - remove cells with no B cell genes (discovered in earlier analysis)
  bcell_min <- 0
  
  # Create histograms of all metrics
  plots <- list()
  plots[[1]] <- ggplot(seurat_obj@meta.data, aes(x=nFeature_RNA)) +
    geom_histogram(color=colour, fill=colour, binwidth=10) + xlim(0,4000) +
    geom_vline(xintercept=nFeature_min) +
    geom_vline(xintercept=nFeature_max) +
    ggtitle(paste("min =", round(nFeature_min) , ", max =", round(nFeature_max)))
  plots[[2]] <- ggplot(seurat_obj@meta.data, aes(x=nCount_RNA)) +
    geom_histogram(color=colour, fill=colour, binwidth=50) + xlim(0,20000) +
    geom_vline(xintercept=nCount_min) +
    geom_vline(xintercept=nCount_max) +
    ggtitle(paste("min =", round(nCount_min) , ", max =", round(nCount_max)))
  plots[[3]] <- ggplot(seurat_obj@meta.data, aes(x=percent.mt)) +
    geom_histogram(color=colour, fill=colour, binwidth=0.25) + xlim(0,100)+
    geom_vline(xintercept=mt_max) +
    ggtitle(paste("max =", round(mt_max,5)))
  plots[[4]] <- ggplot(seurat_obj@meta.data,aes(nCount_RNA, nFeature_RNA, color=percent.mt)) +
    geom_point() + ggtitle(sample) +
    geom_vline(xintercept = nCount_min) + geom_vline(xintercept = nCount_max) +
    geom_hline(yintercept = nFeature_min) + geom_hline(yintercept = nFeature_max) +
    scale_color_gradient(low="grey", high=colour)
  
  pdf(filename)
  print ((plots[[4]] | plots[[1]]) / (plots[[2]] | plots[[3]]))
  dev.off()
  
  # Perform the subset based on threshold values
  temp_obj <- subset(seurat_obj, subset = nFeature_RNA >= nFeature_min & nFeature_RNA < nFeature_max &
                       nCount_RNA >= nCount_min & nCount_RNA < nCount_max &
                       percent.mt < mt_max & percent.Bcell > bcell_min)
  return(temp_obj)
}

# HEALTH QC
for (sample in ls(health_CD19_QC_samples)){
  health_CD19_QC_samples[[sample]] <- QC_sample(health_CD19_QC_samples[[sample]], sample, COLOURS$status[['health']],
                                                file.path(out_dir, paste0(sample, '_QC.pdf')))
}
# LUPUS QC
for (sample in ls(lupus_CD19_QC_samples)){
  lupus_CD19_QC_samples[[sample]] <- QC_sample(lupus_CD19_QC_samples[[sample]], sample, COLOURS$status[['lupus']],
                                               file.path(out_dir, paste0(sample, '_QC.pdf')))
}

# Integrate ------------------------------------
qc_samples <- c(health_CD19_QC_samples, lupus_CD19_QC_samples)
norm_samples <- list()

for (sample in ls(qc_samples)) {
  DefaultAssay(qc_samples[[sample]]) <- 'RNA'
  # Normalise RNA
  norm_samples[[sample]] <- NormalizeData(qc_samples[[sample]], verbose = FALSE)
  # Var features RNA
  norm_samples[[sample]] <- FindVariableFeatures(norm_samples[[sample]], nfeatures = 2000, verbose = FALSE)
  
  # Normalise ADT
  norm_samples[[sample]] <- NormalizeData(norm_samples[[sample]], assay = "ADT", normalization.method = "CLR", verbose = F)
}


num_anchor_genes <- 360

anchors <- FindIntegrationAnchors(object.list = norm_samples, anchor.features = num_anchor_genes)
integrated_samples <- IntegrateData(anchorset = anchors)

integrated_samples <- ScaleData(integrated_samples, verbose = FALSE, assay = 'integrated')
integrated_samples <- ScaleData(integrated_samples, verbose = FALSE, assay = 'ADT')
integrated_samples <- RunPCA(integrated_samples)

dims <- 1:20

integrated_samples <- RunUMAP(integrated_samples, reduction.name='umap', dims=dims)

p <- (VlnPlot(integrated_samples, 'nCount_RNA', group.by = 'Sample', pt.size = 0, cols =   COLOURS$Sample) + NoLegend() |
        VlnPlot(integrated_samples, 'nFeature_RNA', group.by = 'Sample', pt.size = 0, cols = COLOURS$Sample) + NoLegend()) /
  (VlnPlot(integrated_samples, 'percent.mt', group.by = 'Sample', pt.size = 0, cols =   COLOURS$Sample) + NoLegend() |
     VlnPlot(integrated_samples, 'percent.rb', group.by = 'Sample', pt.size = 0, cols =   COLOURS$Sample) + NoLegend())
pdf(file.path(out_dir, 'QC_metrics_by_sample.pdf'))
print(p)
dev.off()

# Checking batch correction
p <- DimPlot(integrated_samples, reduction = 'umap', group.by = 'Sample', cols = COLOURS$Sample) + NoLegend() |
  DimPlot(integrated_samples, reduction = 'umap', group.by = 'Sample', cols = COLOURS$Sample, split.by = 'Sample', ncol = 3)
pdf(file.path(out_dir, 'umap_by_sample.pdf'), width = 14)
print(p)
dev.off()

# Classifying ----
integrated_samples[[MANUAL_CLASSIFICATION]] <- 'Uncertain'
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(2)),] <- CLASSIFICATIONS$TRANSITIONAL
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(0, 1, 4, 5, 6, 7, 8, 13, 14)),] <- CLASSIFICATIONS$NAIVE
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(11)),] <- CLASSIFICATIONS$MZB_2
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(3)),] <- CLASSIFICATIONS$MZB_1
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(18)),] <- CLASSIFICATIONS$DN2
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(9)),] <- CLASSIFICATIONS$ACTIVATED_NAIVE
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(20)),] <- CLASSIFICATIONS$PLASMABLASTS
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(15, 10)),] <- CLASSIFICATIONS$IGM_ONLY
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(12, 16, 17)),] <- CLASSIFICATIONS$CSM
integrated_samples[[MANUAL_CLASSIFICATION]][WhichCells(integrated_samples, idents = c(19)),] <- CLASSIFICATIONS$MALAT1

# Markers ----------------------------------------
integrated_samples@active.ident <- as.factor(integrated_samples$manual_classification)

de_pos_list <- {}
de_neg_list <- {}

for (i in levels(integrated_samples)){
  condition.diffgenes <- FindMarkers(integrated_samples, ident.1 = i)
  pos_dat <- subset(condition.diffgenes, avg_logFC > 0)
  neg_dat <- subset(condition.diffgenes, avg_logFC < 0)
  pos_marker <- rownames(pos_dat)
  neg_marker <- rownames(neg_dat)
  de_pos_list[[i]] <- pos_marker
  de_neg_list[[i]] <- neg_marker
}
de_pos_summary <- t(plyr::ldply(de_pos_list, rbind, .id = "cluster"))
de_neg_summary <- t(plyr::ldply(de_neg_list, rbind, .id = "cluster"))
file_name <- paste0("pos_degenes_summary.csv")
write.csv(de_pos_summary, file=file.path(out_dir, file_name), col.names = F)
file_name <- paste0("neg_degenes_summary.csv")
write.csv(de_neg_summary, file=file.path(out_dir, file_name), col.names = F)

