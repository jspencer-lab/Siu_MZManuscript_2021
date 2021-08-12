#Figure 1 Cydar

#Libraries
library(edgeR)
library(stats)
library(dplyr)
library(viridis)
library(tidyverse)
library(cydar)

#Directories
dir <- file.path("Cydar")
out_dir <- file.path("figure images", "Fig1")

#Import data
count_dat <- readRDS(file.path(dir, "count_dat_Bcellcydar_range.rds"))
y <- DGEList(assay(count_dat), lib.size=count_dat$totals)

#Determining hydrospheres to keep
keep_size <- 100
threshold <- aveLogCPM(keep_size, lib.size=mean(y$samples$lib.size))
keep <- aveLogCPM(y) >= threshold

count_dat_keep <- count_dat[keep,]
y_keep <- y[keep,]
summary(keep)

#Establishing conditions
conditions <- factor(c("Appendix", "mLN", "Spleen", "Appendix", "Spleen", "Appendix", "mLN", "Spleen", "Appendix", "mLN", "Spleen", "Appendix", "mLN", "Spleen", "Appendix", "Spleen"))

design.use <- model.matrix(~ 0 + conditions) 
colnames(design.use) <- levels(conditions)
y_design <- estimateDisp(y_keep, design.use) 

#fit for dispersion using QL
fit.dispersion <- glmQLFit(y_design, design = design.use, robust = TRUE) 

#File names
file_names <- colnames(count_dat)

#### Fig1a
sig.data <- intensities(count_dat_keep)
reranges <- intensityRanges(count_dat_keep, p=0.05)

require(uwot)
umap.out <- umap(sig.data, n_neighbors = 200, min_dist = 0.5)

filename <- "UMAP_Marker.pdf"
png(file.path(outdir, filename), width=20, height=20, units="in", res=300)

lmat <- cbind(matrix(seq_len(6*6), ncol=6, nrow=6), 37)
layout(lmat, widths=c(rep(1, 6), 0.2))
for (i in order(colnames(sig.data))) {
  par(mar=c(2.1, 2.1, 2.1, 2.1))
  out <- plotSphereIntensity(umap.out[,1], umap.out[,2], sig.data[,i], 
                             irange=reranges[,i], main=colnames(sig.data)[i], 
                             xlab="UMAP1", ylab="UMAP2")
}

par(mar=c(0,0,0,0))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(out))
interval <- diff(start.loc)[1]
rect(-0.5, start.loc, 0.5, start.loc+interval, col=out, border=NA)
text(0, -0.5, pos=1, "Low", cex=1.5)
text(0, 0.5+interval,  pos=3, "High", cex=1.5)
text(-0.9, 0, pos=1, srt=90, "Marker intensity", cex=1.5)
dev.off()

#### Fig1b (scaled)
con <- makeContrasts(app_spl = Appendix - Spleen, app_mln = Appendix - mLN, mln_spl = mLN - Spleen, levels = design.use)

res.dispersion <- glmQLFTest(fit.dispersion, contrast = con)

res.dispersion$table$FDR <- spatialFDR(intensities(count_dat_keep), res.dispersion$table$PValue)
a <- 0.05
is.sig <- res.dispersion$table$FDR <= a
summary(is.sig)

sig.coords <- intensities(count_dat_keep)[is.sig,]
sig.res <- res.dispersion[is.sig,]
sig.out <- count_dat_keep[is.sig,]

#Nonredundancy
nonredundant <- findFirstSphere(intensities(sig.out), sig.res$table$PValue)
summary(nonredundant)

important.res <- sig.res[nonredundant,]

all_clusters <- as.numeric(rownames(res.dispersion$table))
cluster_num_sig <- as.numeric(rownames(important.res$table))

cluster_true_false_nonredundant <- all_clusters %in% cluster_num_sig

umap.dat <- data.frame(all_clusters, "UMAP1" = umap.out[,1], "UMAP2" = umap.out[,2], cluster_true_false_nonredundant, is.sig)

counts_all <- y_keep$counts
colnames(counts_all) <- c("292B Appendix", "292B mLN", "292B Spleen", "298C Appendix", "298C Spleen", "318B Appendix", "318B mLN", "318B Spleen", "325C Appendix", "325C mLN", "325C Spleen", "317B Appendix", "317B mLN", "317B Spleen", "343B Appendix", "343B Spleen")

dat_counts <- stack(t(counts_all))
colnames(dat_counts) <- c("file_name", "hydrasphere", "count")

donor <- sub(" .*", "", dat_counts$file_name)
tissue <- sub(".* ", "", dat_counts$file_name)

dat_counts <- data.frame(donor, tissue, dat_counts)

#splitting dat_counts by tissue first then taking hydrasphere means
app_dat.counts <- select(filter(dat_counts, tissue == "Appendix"), c(hydrasphere, count))
spl_dat.counts <- select(filter(dat_counts, tissue == "Spleen"), c(hydrasphere, count))
mln_dat.counts <- select(filter(dat_counts, tissue == "mLN"), c(hydrasphere, count))

#Taking means
app_agg <- app_dat.counts %>% group_by(hydrasphere) %>% summarise(avg=median(count))
spl_agg <- spl_dat.counts %>% group_by(hydrasphere) %>% summarise(avg=median(count))
mln_agg <- mln_dat.counts %>% group_by(hydrasphere) %>% summarise(avg=median(count))

#Combinging (tissues_count is data frame with average count per hydrasphere (row) in y_keep per tissue (column))
tissues_count <- data.frame(app_agg$avg, spl_agg$avg, mln_agg$avg)
colnames(tissues_count) <- c("Appendix_count", "Spleen_count", "mLN_count")

markers <- colnames(sig.data)
tissues_name <- c("Appendix", "Spleen", "mLN")

#Scale
scale_tissues_count <- t(scale(t(tissues_count), center = TRUE))

#Sum of raw
tissues_count$sum <- rowSums(tissues_count)

dat.with.colours <- data.frame(umap.dat, tissues_count)

pdf(file.path(out_dir, paste0("HypersphereTotalCounts_UMAP.pdf")), width=7, height=5)

a <- ggplot(data = dat.with.colours %>% arrange(is.sig, sum), aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(fill = sum, colour = is.sig), size = 1, pch = 21)  + 
  scale_fill_viridis(option = "viridis") + 
  scale_colour_manual(values = c("transparent", "black")) + 
  theme_pubr() +
  labs(x = "UMAP_1", y = "UMAP_2", fill = "Total Count") +
  guides(colour = FALSE) +
  theme(legend.position="right", text = element_text(size = 20), axis.ticks=element_blank(), axis.text=element_blank())

print(a) 

dev.off()


#Plot Scaled by tissue
for (i in 1:3){
  count_use <- scale_tissues_count[,i]
  dat.with.colours <- data.frame(umap.dat, count_use)
  
  pdf(file.path(out_dir, paste0("HypersphereScaledCounts_UMAP_sig05_", tissues_name[i], ".pdf")), width=7, height=5)
  
  a <- ggplot(data = dat.with.colours %>% arrange(is.sig), aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(fill = count_use, size = is.sig, colour = is.sig), pch = 21)  + 
    scale_fill_viridis(option = "cividis") + 
    scale_colour_manual(values = c("transparent", "black")) + 
    theme_pubr() +
    labs(x = "UMAP_1", y = "UMAP_2", fill = "Scaled Count", title = paste0(tissues_name[i]), size = "p < 0.05") +
    guides(colour = FALSE) +
    theme(legend.position="right", text = element_text(size = 20), axis.ticks=element_blank(), axis.text=element_blank()) + 
    scale_size_manual(values = c(0.5, 7)) 
  
  print(a) 
  
  dev.off()
  
}

for (i in 1:3){
  count_use <- scale_tissues_count[,i]
  dat.with.colours <- data.frame(umap.dat, count_use)
  
  pdf(file.path(out_dir, paste0("HypersphereScaledCounts_UMAP_sig05_", tissues_name[i], "_nolegend.pdf")), width=5, height=5.2)
  
  a <- ggplot(data = dat.with.colours %>% arrange(is.sig), aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(fill = count_use, size = is.sig, colour = is.sig), pch = 21)  + 
    scale_fill_viridis(option = "cividis") + 
    scale_colour_manual(values = c("transparent", "black")) + 
    theme_pubr() +
    labs(x = "UMAP_1", y = "UMAP_2", fill = "Scaled Count", title = paste0(tissues_name[i]), size = "p < 0.05") +
    guides(colour = FALSE) +
    theme(legend.position="none", text = element_text(size = 20), axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank(), line = element_blank()) + 
    scale_size_manual(values = c(0.5, 7)) 
  
  print(a) 
  
  dev.off()
  
}

#### Fig1c (heatmap)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dendextend)

dat_phenotype <- data.frame(all_clusters, sig.data)

dat_phenotype_sigdatatf <- data.frame(dat_phenotype, is.sig)
rownames(dat_phenotype_sigdatatf) <- dat_phenotype_sigdatatf$all_clusters


sig_dat_phenotype <- subset(dat_phenotype_sigdatatf, is.sig == TRUE, select = all_clusters:PD.L1)

sig_cluster_id <- sig_dat_phenotype$all_clusters
rownames(sig_dat_phenotype) <- sig_cluster_id

sig_dat_phenotype <- as.matrix(sig_dat_phenotype)

#Take out CD86, PD.L1, FcRL5, FcRL4, CD11c, CD138
take_outheatmap <- c(match(c("CD86", "PD.L1", "FcRL5", "FcRL4", "CD11c", "CD138"), colnames(sig_dat_phenotype)), 
                     grep("^BC-[0-9]", colnames(sig_dat_phenotype)))

takeout_numbers <- setdiff(c(2:27), take_outheatmap)

dat_heatmap <- sig_dat_phenotype[,takeout_numbers]

#Counts
y_keep_sig <- y_keep[is.sig,]

counts_sig <- y_keep_sig$counts
colnames(counts_sig) <- c("292B Appendix", "292B mLN", "292B Spleen", "298C Appendix", "298C Spleen", "318B Appendix", "318B mLN", "318B Spleen", "325C Appendix", "325C mLN", "325C Spleen", "317B Appendix", "317B mLN", "317B Spleen", "343B Appendix", "343B Spleen")

dat_counts <- stack(t(counts_sig))
colnames(dat_counts) <- c("file_name", "hydrasphere", "count")

donor <- sub(" .*", "", dat_counts$file_name)
tissue <- sub(".* ", "", dat_counts$file_name)

dat_counts <- data.frame(donor, tissue, dat_counts)

mean_by_tishy <- dat_counts %>% group_by(hydrasphere, tissue) %>% summarise(average.counts = mean(count))

wide_mean <- mean_by_tishy %>% spread(tissue, average.counts)

ha_dat <- data.frame("hydrasphere" = wide_mean[1], "appendix" = wide_mean[2], "spleen" = wide_mean[4], "mLN"= wide_mean[3])
row.names(ha_dat) <- ha_dat$hydrasphere
ha_dat <- ha_dat[,2:4]


#
dat_heatmap <- data.frame(dat_heatmap, ha_dat)

row_dend = as.dendrogram(hclust(dist(dat_heatmap[,1:20])))
row_dend = rotate(row_dend, 27:32)
row_dend = color_branches(row_dend, k = 6) 


ht <- Heatmap(dat_heatmap[,1:20], col = viridis(100), name = 'Marker\nExpression', cluster_rows = row_dend, show_column_dend = FALSE, row_dend_gp = gpar(lwd = 2), row_dend_width = unit(2, "cm"), column_names_gp = gpar(fontsize = 10), row_split = 6, heatmap_legend_param = list(direction = "horizontal"), column_names_rot = 45, show_row_names = FALSE)

ht


ht_row <- row_order(ht) #[1] 27 26 52  6  8 21 64  9 14 60 67  7 17 29 73 38 33 40 35 44 37 39 41 45 34 36  5 43 75  2  4  3 74 42 46 47 48 11 24 51 12 20 25 18 55 71 28 31 66 54 19 22 10 [54] 32 53  1 72 49 57 23 59 61 58 50 65 69 70 13 68 30 63 15 62 56 16
ht_col <- column_order(ht)

#scaled
scale_dat <- scale(dat_heatmap[,1:20], center = TRUE)

col_fun = colorRamp2(c(-1.5, 0, 1.5), cividis(3))
col_fun(seq(-3, 3))

sht <- Heatmap(scale_dat, col = col_fun, name = 'Column Scaled\nMarker Expression', cluster_rows = row_dend, cluster_columns = FALSE, column_order = ht_col, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 10), row_split = 6, show_row_dend = FALSE, heatmap_legend_param = list(direction = "horizontal"), column_names_rot = 45, border = TRUE, row_names_side = "left")
sht

sht <- Heatmap(scale_dat, col = col_fun, name = 'Column Scaled\nMarker Expression', cluster_rows = row_dend, cluster_columns = FALSE, column_order = ht_col, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 10), row_split = 6, show_row_dend = TRUE, heatmap_legend_param = list(direction = "horizontal"), column_names_rot = 45, border = TRUE, show_row_names = FALSE, row_dend_width = unit(2, "cm"), row_dend_gp = gpar(lwd = 2))

#Plot
pdf(file.path(out_dir, paste0("SigHypersphereHeatmapScaled.pdf")), width=4, height=6)

draw(sht, heatmap_legend_side = "top", annotation_legend_side = "top", auto_adjust = FALSE)

dev.off()

#### Fig 1d (overlay clusterID onto UMAP)
hyperspherename <- rownames(dat_heatmap)
sighypersphere_dat <- {}

for(i in 1:length(ht_row)){
  tempdat <- cbind(ht_row[[i]], rep(i, length(ht_row[[i]])))
  sighypersphere_dat <- rbind(sighypersphere_dat, tempdat)
}

colnames(sighypersphere_dat) <- c("row_id", "subset_id")
sighypersphere_dat <- as.data.frame(sighypersphere_dat)

sighypersphere_dat <- sighypersphere_dat %>% 
                          select(row_id, subset_id) %>% 
                          mutate("all_clusters" = hyperspherename[row_id])
sighypersphere_dat <- sighypersphere_dat[,-1]

for (i in 1:3){
  count_use <- scale_tissues_count[,i]
  dat.with.colours <- data.frame(umap.dat, count_use)
  dat_subset <- plyr::join(dat.with.colours, sighypersphere_dat, by = "all_clusters", type = "left")
  dat_subset$subset_id <- replace_na(dat_subset$subset_id, 0)
  
  #Overlay only (outline clusters only)
  leg_order <- c("1", "2", "3", "4", "5", "6", "0")
  
  pdf(file.path(out_dir, paste0("HypersphereClusterIDlegend_UMAP.pdf")), width=10, height=10)
  
  a <- ggplot(data = dat_subset %>% arrange(is.sig), aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(size = is.sig, fill = as.factor(subset_id), colour = as.factor(subset_id), pch = as.factor(subset_id)), stroke = 1)  + 
    #scale_fill_viridis(option = "cividis") +
    scale_shape_manual(values = c(21, 21, 21, 24, 22, 21, 21)) +
    scale_fill_manual(breaks = leg_order, values = c("0" = "grey", "1" = "#CC476B", "2"=  "#9F7000","3" ="#6CDC52", "4" = "#009681","5" ="#0082CE", "6"= "#B646C7"), labels = c(1, 2, 3, 4, 5, 6, " ")) +
    scale_colour_manual(breaks = leg_order, values = c("0" = "grey", "1" = "black", "2"=  "black","3" ="black", "4" = "black","5" ="black", "6"= "black"), labels = c(1, 2, 3, 4, 5, 6, " ")) + 
    theme_pubr() +
    labs(x = "UMAP_1", y = "UMAP_2", size = "p < 0.05", colour = "Cluster") +
    theme(legend.position="none", text = element_text(size = 25), axis.ticks=element_blank(), axis.text=element_blank()) + 
    scale_size_manual(values = c(0.5, 7)) 
  
  print(a) 
  
  dev.off()
}