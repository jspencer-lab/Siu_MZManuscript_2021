#BCRtisues_supplemental
#For this section, all the analysis was done using the Kleinstein Immcantation (Version 4.0) docker image and RStudio. For more details on how to install this, please visit their documentation: https://immcantation.readthedocs.io/en/stable/docker/intro.html. Majority of this analysis was done prior to the latest Immcantation update which incorporates more streamline 10X analysis into their pipeline.
#https://kleinstein.bitbucket.io/tutorials/intro-lab/index.html

##Libraries
library(plyr) # plyr is needed for circo plot but need to be load before dplyr
library(tidyverse)
library(tidyr)
library(dplyr)
library(tibble)
library(alakazam) #Immcantation
library(shazam) #Immcantation
library(airr) #Immcantation
library(rstatix)
library(kableExtra)
library(glue)
library(circlize)
library(data.table)
library(foreach)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library("ggpubr") #make nicer figures faster
library(patchwork)
library(ggplot2)
library(viridis)
library(randomcoloR) #generate random colours

options(future.globals.maxSize = 4000 * 1024^2)

#Theme
color_use <- c("#7B52DD", "#7DE3D4", "#DC694C", "#77E595", "#E4DB61", "#D78AD6", "#CFDFDE", "#E0E2AC", "#9B9C59", "#778BD5", "#95EA4F", "#D3B6DB", "#75ABC5", "#DA4BA0", "#D49794")

#Directories
out_dir <- file.path("figure images", "Supplemental")

##Set up
db_all_isotype <- subset(db_all, c_call %in% c("IGHM","IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2"))

#Rename isotypes
db_all_isotype$isotype <- forcats::fct_collapse(db_all_isotype$c_call, IGHG = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), IGHA = c("IGHA1", "IGHA2"))


##Abundance
groupby <- "tissue"
a <- estimateAbundance(db_all_isotype, group=groupby, uniform = FALSE)

pdf(file=file.path(out_dir, paste0("ClonalAbundanceCurve_", groupby, ".pdf")),width=3, height=2.5)
plot(a) 
dev.off()

#subset names
groupby <- "subset_names"
a <- estimateAbundance(db_all_isotype, group=groupby, uniform = FALSE, ci =0)

pdf(file=file.path(out_dir, paste0("ClonalAbundanceCurve_noCI_", groupby, ".pdf")),width=3, height=4)
plot(a, col = color_use) #col = color_use for subset names
dev.off()

##Diversity
groupby <- "isotype"
d <- alphaDiversity(db_all_isotype, group=groupby)
p <- plot(d, silent=T)

a <- p + 
  geom_vline(xintercept=c(0,1,2), color="grey50", linetype="dashed") + 
  geom_text(data=data.frame(q=c(0,1,2), y=round(max(p$data$d_upper)/2), 
                            label=c("Richness", "Shannon", "Simpson")), 
            aes(x=q, y=y,label=label), size=3, angle=90, vjust=-0.4, inherit.aes = F, color="grey50")

pdf(file=file.path(out_dir, paste0("HillDiversityCurve_", groupby, ".pdf")),width=5, height=10)
plot(a)
dev.off()

#subset
groupby <- "subset_names"
d <- alphaDiversity(db_all_isotype, group=groupby, ci = 0)
p <- plot(d, silent=T, col = color_use)

a <- p + 
  geom_vline(xintercept=c(0,1,2), color="grey50", linetype="dashed") + 
  geom_text(data=data.frame(q=c(0,1,2), y=round(max(p$data$d_upper)/1.15), 
                            label=c("Richness", "Shannon", "Simpson")), 
            aes(x=q, y=y,label=label), size=3, angle=90, vjust=-0.4, inherit.aes = F, color="grey50")

pdf(file=file.path(out_dir, paste0("HillDiversityCurve_noCI_", groupby, ".pdf")),width=5, height=5)
plot(a)
dev.off()

## Calculate CDR3 amino acid properties
db_aa <- aminoAcidProperties(db_all_isotype, seq="junction", nt=T, trim=T, label="cdr3")
# Plot
groupby <- "isotype"
a <- ggplot(db_aa, aes(x=isotype, y=cdr3_aa_length)) + theme_bw() +
  ggtitle("CDR3 length") + xlab("Isotype") + ylab("Amino acids") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=isotype))

pdf(file=file.path(out_dir, paste0("CDR3length_", groupby, ".pdf")),width=6, height=5)
plot(a)
dev.off()

groupby <- "tissue"
a <- ggplot(db_aa, aes(x=tissue, y=cdr3_aa_length)) + theme_bw() +
  ggtitle("CDR3 length") + xlab("Tissue") + ylab("Amino acids") +
  scale_fill_manual(name="Tissue", values=c("blue", "black", "red")) +
  geom_boxplot(aes(fill=tissue))

pdf(file=file.path(out_dir, paste0("CDR3length_", groupby, ".pdf")),width=6, height=5)
plot(a)
dev.off()

groupby <- "subset_names"
a <- ggplot(db_aa, aes(x=subset_names, y=cdr3_aa_length)) + theme_bw() +
  ggtitle("CDR3 length") + xlab("Subset names") + ylab("Amino acids") +
  scale_fill_manual(name="Subset names", values=color_use) +
  geom_boxplot(aes(fill=subset_names))

pdf(file=file.path(out_dir, paste0("CDR3length_", groupby, ".pdf")),width=15, height=5)
plot(a)
dev.off()

##SHM targetting bias
# Build and plot SHM targeting model
m <- createTargetingModel(db_aa, vCallColumn="v_call")
# nucleotides: center nucleotide characters to plot
plotMutability(m, nucleotides=c("A","C"), size=1.2)

## TS and Naive clones (isolate cells that are naive or TS)

db_tn <- subset(db_all, subset_names %in% c("Naive", "TS"))

dat_list <- db_tn %>% group_by(donor, clone_id) %>% summarise(n = n()) %>% ungroup()
dat_list <- subset(dat_list, n > 1)

cell_id_list <- data.frame("subset_names" = db_tn$subset_names, "cell_id" = db_tn$cell_id, "donor" = db_tn$donor, "clone_id" = db_tn$clone_id)
dat_list <- merge(dat_list, cell_id_list, by = c("donor", "clone_id"))
dat_list <- dat_list[order(-dat_list$n),]

write.csv(dat_list, file = file.path(out_dir, "TS_naiveclones.csv"))

#Fig S3H MZB clone UMAP
donor_list <- c("A", "B", "C")
db_all <- {}

for(i in 1:3){
  donor <- donor_list[i]
  donor <- "C"
  db <- readChangeoDb(file.path(in_dir, paste0("donor", donor, "_data_ph_genotyped_germ-pass.tsv")))
  
  #Import seurat_metadata
  #seurat_metadata <- readRDS(paste0("C:/Users/jacqu/Documents/JS10X_013_ALL/data/integratedall_metadata_donor", donor, ".rds"))
  seurat_metadata <- readRDS(paste0("C:/Users/Jacqueline Siu/OneDrive - King's College London/Science/Projects/2020_tissues10X/Data/Seurat saved objects/JS10X_013_ALL integratedall_metadata_donor", donor, ".rds"))
  seurat_metadata  <- rownames_to_column(seurat_metadata , "cell_id")
  seurat_metadata <- seurat_metadata[!duplicated(seurat_metadata$cell_id), ]
  
  db <- inner_join(seurat_metadata, db,  by='cell_id')
  
  subset_ids <- data.frame(seurat_clusters = 0:22, 
                           subset_names = c("MZB-1", "TS", "ABC3", "ABC1", "Memory", "aNAV", "Naive", "Memory", "ABC2", "IgM-only", "MZB-2", "Naive", "Memory", "GC", "IgM-only", "GC", "Naive", "GC", "DN-A", "aNAV", "DN-B", "ABC4", "PB"))
  
  subset_ids$seurat_clusters <- as.factor(subset_ids$seurat_clusters)
  
  #Add new column to data
  db <- dplyr::right_join(db, subset_ids, by = c("seurat_clusters"))
  
  db$sequence_id <- db$cell_id
  
  # Calculate combined R and S mutation frequencies
  db <- observedMutations(db, sequenceColumn="sequence_alignment",
                          germlineColumn="germline_alignment_d_mask",
                          regionDefinition=NULL,
                          frequency=T, 
                          combine=T,
                          nproc=16)
  
  db$seurat_clusters <- as.factor(db$seurat_clusters)
  seurat_levels <- c(18, 2, 20, 3, 17, 13, 15, 14, 7, 12, 0, 4, 1, 6, 16, 5, 9, 11, 8, 19, 10, 21, 22)
  db <- transform(db, seurat_clusters = factor(seurat_clusters, levels = seurat_levels))
  
  db$tissue <- as.factor(db$tissue.x)
  
  db <- db %>% dplyr::mutate(mutation = 1 - v_identity)
  
  #For db_all
  #db_all <- db
  db_all <- Map(c, db_all, db) 
}
### END OF LOOP ###

db_all$seurat_clusters <- as.factor(db_all$seurat_clusters)
seurat_levels <- c(18, 2, 20, 3, 17, 13, 15, 14, 7, 12, 0, 4, 1, 6, 16, 5, 9, 11, 8, 19, 10, 21, 22)
db_all <- transform(db_all, seurat_clusters = factor(seurat_clusters, levels = seurat_levels))

db_all$tissue <- as.factor(db_all$tissue.x)

#saved Seurat object
dat <- readRDS(paste0("integrated_tissue.rds"))

umap_dat <- as.data.frame(Embeddings(dat[["umap"]]))

umap_dat  <- rownames_to_column(umap_dat , "cell_id")

db_all$donorclone_id <-paste0(db_all$donor, "_", db_all$clone_id) 

dt <- data.table(db_all)
groupByFormatted <- 'donorclone_id'
dt <- dt[ , list( yidx = list(.I) ) , by=groupByFormatted ] #yidx = indices of repeats
groups <- dt[,yidx] 
groupsLength <- unlist(lapply(groups,length))
groupsSingletons <- groups[groupsLength==1]
groups <- groups[groupsLength>1]
numbOfTotalGroups <- length(groups)

## Determine what clones are disseminate or not
listResults  <- #List of matching clones between subsets / tissues
  foreach(i=iterators::icount(numbOfTotalGroups)) %do% {
    dfClone <- db_all[groups[[i]],c("donor", "clone_id", "tissue", "subset_names", "mu_freq", "UMAP_1", "UMAP_2", "donorclone_id")]
    
    #Resident vs Disseminate
    uniqueCategories <- length( unique( as.character(dfClone$tissue)) )
    if( uniqueCategories==1 ){
      clone_type <- "Resident"
    }else{
      clone_type <- "Disseminate"
    }
    
    #SHM vs no
    uniquemufreq <- sd(dfClone$mu_freq)
    
    if(uniquemufreq==0){
      SHM_type <- "noChangeMutation"
      
    }else{
      SHM_type <- "Mutation"
    }
    
    length_dfClone <- nrow(dfClone)
    if (length_dfClone == 1){
      matReturn <- NULL
    } else if ("MZB-1" %in% dfClone$subset_names & "MZB-2" %in% dfClone$subset_names){
      mzb_type <- "both"
    } else if ("MZB-1" %in% dfClone$subset_names){
      mzb_type <- "MZB1"
    } else if ("MZB-2" %in% dfClone$subset_names){
      mzb_type <- "MZB2"
    } else {
      mzb_type <- "none"
    }
    
    matReturn <- cbind(dfClone, clone_type, SHM_type, mzb_type)
    
    
    return(matReturn)
  } 

df <- do.call(rbind, listResults)
df <- as.data.frame(df)

df_mz <- subset(df, mzb_type != "none" )

cellcounts_dat <- df_mz %>% group_by(donor, mzb_type) %>% summarise(cellcount_allsubsets = n())

db_all <- dplyr::left_join(db_all, umap_dat, by = c("cell_id"))

df_mz <- subset(df, subset_names == "MZB-1" | subset_names == "MZB-2" )

plot <- ggplot(df_mz, aes(x = UMAP_1, y = UMAP_2, colour = mzb_type, shape = subset_names)) +
  geom_point() +
  #facet_grid(tissue ~ donor) +
  theme_classic2() +
  labs(colour = "Clone contains...") +
  scale_colour_manual(values = c("orange3", "#95EA4F", "#D3B6DB")) +
  theme(legend.position = "none")

file_name <- paste0("MZBclonaloverlay_all_combined_nolegend.png")
ggsave(filename = file.path(out_dir, file_name), plot = plot, width = 10, height = 10, units = c("cm"))
