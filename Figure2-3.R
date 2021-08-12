#Figure2-3 (Figure2d, Figure 3a, c, d)
#scVDJ analysis of tisues

#For this section, all the analysis was done using the Kleinstein Immcantation (Version 4.0) docker image and RStudio. For more details on how to install this, please visit their documentation: https://immcantation.readthedocs.io/en/stable/docker/intro.html. Majority of this analysis was done prior to the latest Immcantation update which incorporates more streamline 10X analysis into their pipeline.

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
out_dir <- file.path("figure images", "Fig3")

in_dir <- file.path("data", "results", "changeo")

seurat_metadata <- readRDS(paste0("/Data/Seurat saved objects/integratedall_metadata_donorB.rds"))


#Importing, integrating VDJ-Seurat, pooling donor data
donor_list <- c("A", "B", "C")
db_all <- {}

for(i in 1:3){
donor <- donor_list[i]
donor <- "C"
db <- readChangeoDb(file.path(in_dir, paste0("donor", donor, "_data_ph_genotyped_germ-pass.tsv")))

#Import seurat_metadata
seurat_metadata <- readRDS(paste0("data/integratedall_metadata_donor", donor, ".rds"))
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

#### Fig 2d
db_all <- db_all %>% dplyr::mutate(mutation = 1 - v_identity)

g1 <- db_all  %>%
  select(subset_names, tissue, mu_freq) %>% # !!as.name(i) allow to pass a string sotred in a variable as a column name
  group_by(subset_names) %>%
  ggplot( aes(x=subset_names, y=mu_freq, fill = subset_names))+
  geom_boxplot(notch = FALSE, varwidth = FALSE) +
  scale_fill_manual(name="Subsets", values=color_use, guide=FALSE) +
  #facet_grid(cols = vars(tissue)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), text = element_text(size=20), axis.text.x=element_text(angle=45,hjust=1), legend.position = "none")+
  ylab("Mutation Frequency") +
  ylim(0, 0.16)
ggsave(plot = g1, filename = paste0(out_dir, "/MutationFrequency_allDonors.pdf"), width = 20, height = 15, unit = c("cm"))

#Stats
res.aov <- aov(mu_freq ~ subset_names, data = db_all) #Levenes signigifcant, not normal
summary(res.aov)
TukeyHSD(res.aov)

plot(res.aov, 1)
library(car)
leveneTest(mu_freq ~ subset_names, data = db_all)
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov )
shapiro.test(x = aov_residuals)

kruskal.test(mu_freq ~ subset_names, data = db_all)


#### Fig 3a (dissemination ratio)
statReturn <- NULL
donor_list <- c("A", "B", "C")

## Import Data (do the big loop here) 
for (i in 1:3){
  donor <- donor_list[i]
  
  #import 
  db <- readChangeoDb(file.path(in_dir, paste0("donor", donor, "_data_ph_genotyped_germ-pass.tsv")))
  
  #Import seurat_metadata
  seurat_metadata <- readRDS(paste0("integratedall_metadata_donor", donor, ".rds"))
  seurat_metadata  <- rownames_to_column(seurat_metadata , "cell_id")
  seurat_metadata <- seurat_metadata[!duplicated(seurat_metadata$cell_id), ]
  
  #join 
  db_join <- inner_join(seurat_metadata, db,  by='cell_id')
  db_join$sequence_id <- db_join$cell_id
  
  #make dataframe for circos
  db <- db_join
  db <- as.data.frame(db)
  
  ##Adding simplified subset names
  subset_ids <- data.frame(seurat_clusters = 0:22, 
                           subset_names = c("MZ-A", "TS", "MALAT1", "ActB1", "CD27+IgG", "ANAV", "Naive", "CD27+IgA", "ActB2", "IgM-only", "MZ-B", "Naive", "CD27+IgG", "GC", "IgM-only", "GC", "Naive", "GC", "DN-A", "ANAV", "DN-B", "Interferon", "PB"))
  
  subset_ids$seurat_clusters <- as.factor(subset_ids$seurat_clusters)
  
  #Add new column to data
  db <- dplyr::right_join(db, subset_ids, by = c("seurat_clusters"))
  
  # Determine clones with 2+ cells that are tissue-resident vs disseminate
  db$tissue <- factor(db$tissue.x, levels=c("APP","MLN","SPL"), labels=c("APP","MLN","SPL"))
  db$seurat_clusters_old <- factor(db$seurat_clusters)
  
  db$seurat_clusters <- factor(db$subset_names)
  
  db <- observedMutations(db, sequenceColumn="sequence_alignment",
                          germlineColumn="germline_alignment_d_mask",
                          regionDefinition=NULL,
                          frequency=T, 
                          combine=T,
                          nproc=16)
  
  dt <- data.table(db)
  groupByFormatted <- 'clone_id'
  dt <- dt[ , list( yidx = list(.I) ) , by=groupByFormatted ] #yidx = indices of repeats
  groups <- dt[,yidx] 
  groupsLength <- unlist(lapply(groups,length))
  groupsSingletons <- groups[groupsLength==1]
  groups <- groups[groupsLength>1]
  numbOfTotalGroups <- length(groups)
  
  ## Determine what clones are disseiminate or not
  listResults  <- #List of matching clones between subsets / tissues
    foreach(i=iterators::icount(numbOfTotalGroups)) %do% {
      dfClone <- db[groups[[i]],c("clone_id", "tissue", "subset_names",
                                  "mu_freq")]
      cloneID <- dfClone$clone_id[1]
      
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
      
      rows <- nrow(dfClone)
      matReturn <- NULL
      
      for(a in 1:rows){
        tissue_id <- dfClone[a,"tissue"]
        subset_id <- dfClone[a,"subset_names"]
        
        
        matReturn = rbind(matReturn, c("cloneID"= cloneID,
                                       "tissue"= tissue_id,
                                       "subset" = subset_id,
                                       "clone_type" = clone_type, 
                                       "SHM_type" = SHM_type))
        
        
      }
      return(matReturn)
    } 
  
  df <- do.call(rbind, listResults)
  df <- as.data.frame(df)
  
  #Make summary stats
  subset_names <- unique(db$subset_names)
  
  for (i in 1:length(subset_names)){
    subset_use <- subset_names[i]
    temp_dat <- subset(df, subset == subset_use, select = c("cloneID", "clone_type", "SHM_type"))
    temp_dat_unique <- unique(temp_dat)
    
    temp_residentcount <- length(which(temp_dat_unique$clone_type == "Resident"))
    temp_disscount <- length(which(temp_dat_unique$clone_type == "Disseminate"))
    temp_percent <- (temp_disscount/temp_residentcount)*100
    
    statReturn <- rbind(statReturn, c("Donor" = donor, 
                                      "Subset" = subset_use,
                                      "Count_Res" = as.numeric(temp_residentcount),
                                      "Count_Dis" = as.numeric(temp_disscount),
                                      "Percent" = as.numeric(temp_percent)))
  }
}

statReturn <- NULL
donor_list <- c("A", "B", "C")

## Import Data (do the big loop here) 
for (i in 1:3){
  donor <- donor_list[i]
  
  #import 
  db <- readChangeoDb(paste0("data/results/changeo/donor", donor, "_data_ph_genotyped_germ-pass.tsv"))
  
  #Import seurat_metadata
  seurat_metadata <- readRDS(paste0("data/integratedall_metadata_donor", donor, ".rds"))
  seurat_metadata  <- rownames_to_column(seurat_metadata , "cell_id")
  seurat_metadata <- seurat_metadata[!duplicated(seurat_metadata$cell_id), ]
  
  #join 
  db_join <- inner_join(seurat_metadata, db,  by='cell_id')
  db_join$sequence_id <- db_join$cell_id
  
  #make dataframe for circos
  db <- db_join
  db <- as.data.frame(db)
  
  ##Adding simplified subset names
  subset_ids <- data.frame(seurat_clusters = 0:22, 
                           subset_names = c("MZB-1", "TS", "ABC3", "ABC1", "Memory", "aNAV", "Naive", "Memory", "ABC2", "IgM-only", "MZB-2", "Naive", "Memory", "GC", "IgM-only", "GC", "Naive", "GC", "DN-A", "aNAV", "DN-B", "ABC4", "PB"))
  
  subset_ids$seurat_clusters <- as.factor(subset_ids$seurat_clusters)
  
  #Add new column to data
  db <- dplyr::right_join(db, subset_ids, by = c("seurat_clusters"))
  
  # Determine clones with 2+ cells that are tissue-resident vs disseminate
  db$tissue <- factor(db$tissue.x, levels=c("APP","MLN","SPL"), labels=c("APP","MLN","SPL"))
  db$seurat_clusters_old <- factor(db$seurat_clusters)
  
  db$seurat_clusters <- factor(db$subset_names)
  
  db <- observedMutations(db, sequenceColumn="sequence_alignment",
                          germlineColumn="germline_alignment_d_mask",
                          regionDefinition=NULL,
                          frequency=T, 
                          combine=T,
                          nproc=16)
  
  dt <- data.table(db)
  groupByFormatted <- 'clone_id'
  dt <- dt[ , list( yidx = list(.I) ) , by=groupByFormatted ] #yidx = indices of repeats
  groups <- dt[,yidx] 
  groupsLength <- unlist(lapply(groups,length))
  groupsSingletons <- groups[groupsLength==1]
  groups <- groups[groupsLength>1]
  numbOfTotalGroups <- length(groups)
  
  ## Determine what clones are disseiminate or not
  listResults  <- #List of matching clones between subsets / tissues
    foreach(i=iterators::icount(numbOfTotalGroups)) %do% {
      dfClone <- db[groups[[i]],c("clone_id", "tissue", "subset_names",
                                  "mu_freq")]
      cloneID <- dfClone$clone_id[1]
      
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
      
      rows <- nrow(dfClone)
      matReturn <- NULL
      
      for(a in 1:rows){
        tissue_id <- dfClone[a,"tissue"]
        subset_id <- dfClone[a,"subset_names"]
        
        
        matReturn = rbind(matReturn, c("cloneID"= cloneID,
                                       "tissue"= tissue_id,
                                       "subset" = subset_id,
                                       "clone_type" = clone_type, 
                                       "SHM_type" = SHM_type))
        
        
      }
      return(matReturn)
    } 
  
  df <- do.call(rbind, listResults)
  df <- as.data.frame(df)
  
  #Make summary stats
  subset_names <- unique(db$subset_names)
  
  for (i in 1:length(subset_names)){
    subset_use <- subset_names[i]
    temp_dat <- subset(df, subset == subset_use, select = c("cloneID", "clone_type", "SHM_type"))
    temp_dat_unique <- unique(temp_dat)
    
    temp_residentcount <- length(which(temp_dat_unique$clone_type == "Resident"))
    temp_disscount <- length(which(temp_dat_unique$clone_type == "Disseminate"))
    temp_percent <- (temp_disscount/temp_residentcount)*100
    
    statReturn <- rbind(statReturn, c("Donor" = donor, 
                                      "Subset" = subset_use,
                                      "Count_Res" = as.numeric(temp_residentcount),
                                      "Count_Dis" = as.numeric(temp_disscount),
                                      "Percent" = as.numeric(temp_percent)))
  }
}

#Normalise between donors
statReturn3 <- as.data.frame(statReturn)
statReturn3$Percent <- as.numeric(statReturn3$Percent)
avgratio_donor <- statReturn3 %>% group_by(Donor) %>% summarise(avg_percent = mean(Percent, na.rm = T))

statReturn3 <- right_join(statReturn3, avgratio_donor)

#Log Change (fold change of Dis/Res ratio over average ratios per donor)
statReturn4 <- statReturn3 %>% mutate(LogChange = Percent/ avg_percent)


stat_plot <- ggplot(statReturn4, aes(x = Subset, y = LogChange)) + 
  geom_boxplot() +
  geom_point(aes(colour = Donor), size = 5) +
  labs(y = "Normalized Dissiminated / Resident Clones") +
  theme_pubr() +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_hline(yintercept=1, linetype="dashed", color = "black")


file_name <- paste0("percentage_disseminatedclones_foldchange.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = stat_plot, width = 25, height = 15, units = c("cm"))

#### Fig 3c (correlation tables)
library(RcmdrMisc) #rcorr.adjust function
library(corrplot) #to make correlation plots

clone_tables <- list()
donor_list <- c("A", "B", "C")
tissue_list <- c("APP", "MLN", "SPL")

for (j in 1:length(donor_list)){
  donor <- donor_list[j]
  db <- readChangeoDb(paste0("data/results/changeo/donor", donor, "_data_ph_genotyped_germ-pass.tsv"))
  
  #Import seurat_metadata
  seurat_metadata <- readRDS(paste0("integratedall_metadata_donor", donor, ".rds"))
  seurat_metadata  <- rownames_to_column(seurat_metadata , "cell_id")
  seurat_metadata <- seurat_metadata[!duplicated(seurat_metadata$cell_id), ]
  
  #join 
  db_join <- inner_join(seurat_metadata, db,  by='cell_id')
  db_join$sequence_id <- db_join$cell_id
  
  #make dataframe for circos
  db <- db_join
  db <- as.data.frame(db)
  
  ##Adding simplified subset names
  #Aug 26 2020 names
  subset_ids <- data.frame(seurat_clusters = 0:22, 
                           subset_names = c("MZB-1", "TS", "ABC3", "ABC1", "Memory", "aNAV", "Naive", "Memory", "ABC2", "IgM-only", "MZB-2", "Naive", "Memory", "GC", "IgM-only", "GC", "Naive", "GC", "DN-A", "aNAV", "DN-B", "ABC4", "PB"))
  
  
  subset_ids$seurat_clusters <- as.factor(subset_ids$seurat_clusters)
  
  #Add new column to data
  db <- dplyr::right_join(db, subset_ids, by = c("seurat_clusters"))
  
  for (i in 1:length(tissue_list)){
    tissue <- tissue_list[i]
    
    g <- db
    g <- g[g$tissue.x == tissue,]
    g <- g[g$clone_id %in% names(which(table(g$clone_id) > 1)),]
    g <- g[,c('clone_id', 'subset_names')]
    
    sample_name <- paste0(donor, "_", tissue)
    
    clone_tables[[sample_name]] <- as.data.frame.matrix(table(g))
    rownames(clone_tables[[sample_name]]) <- paste0(sample_name, '_',
                                                    rownames(clone_tables[[sample_name]]))
    # to change it to T/F if there is a cell (not # of cells per subset)
    for (c in colnames(clone_tables[[sample_name]])){
      clone_tables[[sample_name]][[c]] <- as.integer(clone_tables[[sample_name]][[c]] > 0)
    }
  }
}

#Corr plots for each tissue (combined the donors)
for (x in 1:length(tissue_list)){
  tissue2 <- tissue_list[x]
  tab_tissue <- c(0+x, 3 +x, 6 +x)
  
  tissue_all <- plyr::rbind.fill(clone_tables[[tab_tissue[1]]], clone_tables[[tab_tissue[2]]],
                                 clone_tables[[tab_tissue[3]]])
  tissue_cor <- rcorr.adjust(as.matrix(tissue_all), type="spearman")
  
  #replace all < with 0
  tissue_cor$P[tissue_cor$P == "<.0001"] <- 0.0001
  tissue_cor$P<- data.matrix(tissue_cor$P)
  class(tissue_cor$P) <- "numeric"
  
  classif_order <- c('ABC1', 'ABC2', 'ABC3', 'ABC4', 'aNAV', 'DN-A', 'DN-B', 'GC', 'IgM-only', 'Memory', 'MZB-1', 'MZB-2', 'Naive', 'PB', 'TS')
  
  tissue_cor$R$r <- tissue_cor$R$r[classif_order, classif_order]
  tissue_cor$P <- tissue_cor$P[classif_order, classif_order]
  tissue_cor$R$P <- tissue_cor$R$P[classif_order, classif_order]
  
  
  file_name <- paste0("corr_alphaorder_3siglevels_", tissue2, ".pdf")
  pdf(file = file.path(out_dir, file_name), width = 5, height = 5)
  par(xpd=TRUE)
  corrplot(tissue_cor$R$r, type="lower", order="original", method = "color",
           p.mat = tissue_cor$P, sig.level = c(.001, .01, .05), pch.cex = 1,
           insig = "label_sig", tl.col = "black", tl.srt = 90, outline = "gray",
           title = paste0(tissue2), mar=c(0,0,2,0), diag = FALSE)

  dev.off()
  
  file_name <- paste0("corr_alphaorder_blank_", tissue2, ".pdf")
  pdf(file = file.path(out_dir, file_name), width = 5, height = 5)
  par(xpd=TRUE)
  corrplot(tissue_cor$R$r, type="lower", order="original", method = "color", 
           p.mat = tissue_cor$P, sig.level = c(.001, .01, .05), pch.cex = 1, 
           insig = "n", tl.col = "black", tl.srt = 90, outline = "gray",
           title = paste0(tissue2), mar=c(0,0,2,0), diag = FALSE)
  
  dev.off()
  
}

#### Fig 3d (circos)
donor <- "B"

#import 
db <- readChangeoDb(paste0("data/results/changeo/donor", donor, "_data_ph_genotyped_germ-pass.tsv"))

#Import seurat_metadata
seurat_metadata <- readRDS(paste0("integratedall_metadata_donor", donor, ".rds"))
seurat_metadata  <- rownames_to_column(seurat_metadata , "cell_id")
seurat_metadata <- seurat_metadata[!duplicated(seurat_metadata$cell_id), ]

#join 
db_join <- inner_join(seurat_metadata, db,  by='cell_id')
db_join$sequence_id <- db_join$cell_id

#make dataframe for circos
db <- db_join
db <- as.data.frame(db)

##Adding simplified subset names
subset_ids <- data.frame(seurat_clusters = 0:22, 
                         subset_names = c("MZB-1", "TS", "ABC3", "ABC1", "Memory", "aNAV", "Naive", "Memory", "ABC2", "IgM-only", "MZB-2", "Naive", "Memory", "GC", "IgM-only", "GC", "Naive", "GC", "DN-A", "aNAV", "DN-B", "ABC4", "PB"))


subset_ids$seurat_clusters <- as.factor(subset_ids$seurat_clusters)

#Add new column to data
db <- dplyr::right_join(db, subset_ids, by = c("seurat_clusters"))

# Determine clones with 2+ cells that are tissue-resident vs disseminate
db$tissue <- factor(db$tissue.x, levels=c("APP","MLN","SPL"), labels=c("APP","MLN","SPL"))
db$seurat_clusters_old <- factor(db$seurat_clusters)

db$seurat_clusters <- factor(db$subset_names)

# Create new combination
db$TIS_SC <- paste0(db$tissue, "|", db$seurat_clusters)

# Order
NAMES <- unique(db$TIS_SC)
list_NAMES <- lapply(NAMES, function(x){ strsplit(x, "\\|")[[1]] })
df_NAMES <- do.call(rbind.data.frame, list_NAMES)
colnames(df_NAMES) <- c("tissue", "seurat_clusters")
df_NAMES$tissue <- factor(df_NAMES$tissue , levels=levels(db$tissue))
df_NAMES$seurat_clusters <- factor(df_NAMES$seurat_clusters, levels=levels(db$seurat_clusters))
NAMES_ORDER <- NAMES[order(df_NAMES$tissue, df_NAMES$seurat_clusters)]
db$TIS_SC <- factor(db$TIS_SC, levels=NAMES_ORDER)

#COLOURS
COLORS <- rep(color_use, 3)
nCOLORS <- length(COLORS)

COLORS_NAMES <- levels(db$TIS_SC)
names(COLORS) <- COLORS_NAMES

# Circos plot with unique links (collapsed clones)

#Make data frame of only connections (no value)
table(db$TIS_SC, db$seurat_clusters)
dt <- data.table(db)
groupByFormatted <- 'clone_id'
dt <- dt[ , list( yidx = list(.I) ) , by=groupByFormatted ] #yidx = indices of repeats
groups <- dt[,yidx] 
groupsLength <- unlist(lapply(groups,length))
groups <- groups[groupsLength>1]
numbOfTotalGroups <- length(groups)

listResults2  <- #List of matching clones between subsets / tissues
  foreach( i=iterators::icount(numbOfTotalGroups)) %do% {
    dfClone <- db[groups[[i]],c("clone_id", "TIS_SC", "tissue",
                                "seurat_clusters",
                                "seurat_clusters_old")]
    cloneID <- dfClone$clone_id[1]
    dfClone <- unique(dfClone)
    uniqueTissues <- length(unique(as.character(dfClone$tissue)))
    
    if( uniqueTissues==1 ){
      clone_type <- "Resident"
      
    }else{
      clone_type <- "Disseminate"
      
    }
    
    length_dfClone <- nrow(dfClone)
    
    if (length_dfClone == 1){
      matReturn <- NULL
    } else {
      
      matReturn <- NULL
      rows <- nrow(dfClone)
      for( a in 1:(rows-1) ){
        for( b in (a+1):rows) {
          a_num <- as.integer(dfClone[a, "seurat_clusters_old"])
          b_num <- as.integer(dfClone[b, "seurat_clusters_old"])
          
          a_tis <- as.factor(dfClone[a, "tissue"])
          b_tis <- as.factor(dfClone[b, "tissue"])
          
          
          if (a_num > b_num) {
            from <- dfClone[b,"TIS_SC"]
            to <- dfClone[a,"TIS_SC"]
            pair_name <- paste0(dfClone[b,"TIS_SC"], "_", 
                                dfClone[a,"TIS_SC"])
          } else if(a_num < b_num){ 
            from <- dfClone[a,"TIS_SC"]
            to <- dfClone[b,"TIS_SC"]
            pair_name <- paste0(dfClone[a,"TIS_SC"], "_", 
                                dfClone[b,"TIS_SC"])
          } else if (a_num == b_num & a_tis == "APP"){
            from <- dfClone[a,"TIS_SC"]
            to <- dfClone[b,"TIS_SC"]
            pair_name <- paste0(dfClone[a,"TIS_SC"], "_", 
                                dfClone[b,"TIS_SC"])
          } else if (a_num == b_num & b_tis == "APP") {
            from <- dfClone[b,"TIS_SC"]
            to <- dfClone[a,"TIS_SC"]
            pair_name <- paste0(dfClone[b,"TIS_SC"], "_", 
                                dfClone[a,"TIS_SC"])
          } else if (a_num == b_num & a_tis == "MLN"){
            from <- dfClone[a,"TIS_SC"]
            to <- dfClone[b,"TIS_SC"]
            pair_name <- paste0(dfClone[a,"TIS_SC"], "_", 
                                dfClone[b,"TIS_SC"])
          } else if (a_num == b_num & b_tis == "MLN"){
            from <- as.character(dfClone[b,"TIS_SC"])
            to <- dfClone[a,"TIS_SC"]
            pair_name <- paste0(dfClone[b,"TIS_SC"], "_", 
                                dfClone[a,"TIS_SC"])
          }
          
          matReturn <- rbind(matReturn, c("from" = as.character(from),
                                          "to"= as.character(to), "value" = 1, 
                                          "clone" = cloneID,
                                          "clone_type" = clone_type,
                                          "pair_name" = paste0(pair_name, 
                                                               "_", clone_type))
          )
          
        }
      }
    }
    return(matReturn)
  } 

df_links <- do.call(rbind, listResults2)
df_links <- as.data.frame(df_links) 

#
colnames(df_links) <- c("from", "to", "value", "clone_id", "clone_type", "pair_name")
df_links$value <- as.integer(df_links$value)
df_links$from <- as.factor(df_links$from)
df_links$to <- as.factor(df_links$to)
df_links$clone_type <- as.factor(df_links$clone_type)
df_links$pair_name <- as.factor(df_links$pair_name)

#Disseminate only
df_links <- subset(df_links, df_links$clone_type == "Disseminate")

count <- df_links %>% dplyr::count(pair_name)
count <- as.data.frame(count)
colnames(count) <- c("pair_name", "freq")
count <- count[order(count$freq, decreasing= TRUE),]

#Add back columns to count_freq data
label_dat <- lapply(as.character(count$pair_name), function(x){ strsplit(x, "\\_")[[1]] })
df_label <- do.call(rbind.data.frame, label_dat)
colnames(df_label) <- c("from", "to", "clone_type")
label_clone_type <- df_label$clone_type

label_datfrom <- lapply(df_label$from, function(x){ strsplit(x, "\\|")[[1]] })
df_labelfrom <- do.call(rbind.data.frame, label_datfrom)
colnames(df_labelfrom) <- c("from_tissue", "from_cluster")

label_datto <- lapply(df_label$to, function(x){ strsplit(x, "\\|")[[1]] })
df_labelto <- do.call(rbind.data.frame, label_datto)
colnames(df_labelto) <- c("to_tissue", "to_cluster")

# 2 Circcos Plot
dat_circcos <- data.frame("from" = df_label$from, "to" = df_label$to, "freq" = count$freq, "clone_type" = df_label$clone_type)

dat_circcos$freq <- as.integer(dat_circcos$freq)
dat_circcos$from <- as.factor(dat_circcos$from)
dat_circcos$to <- as.factor(dat_circcos$to)

TIS_SC_COMBO <- unique( c(levels(dat_circcos$from), levels(dat_circcos$to)) )
GRID.COL <- COLORS[ match( TIS_SC_COMBO, names(COLORS) ) ]

list_NAMES <- lapply(TIS_SC_COMBO, function(x){ strsplit(x, "\\|")[[1]] })
df_NAMES <- do.call(rbind.data.frame, list_NAMES)
colnames(df_NAMES) <- c("TISSUE", "SEURAT_CLUSTERS")
df_NAMES$TISSUE <- factor(df_NAMES$TISSUE )
df_NAMES$SEURAT_CLUSTERS <- as.factor(df_NAMES$SEURAT_CLUSTERS)
df_NAMES$SEURAT_CLUSTERS <- factor(df_NAMES$SEURAT_CLUSTERS,
                                   levels=levels(df_NAMES$SEURAT_CLUSTERS))

df_NAMES$SEURAT_CLUSTERS <- factor(df_NAMES$SEURAT_CLUSTERS,
                                   levels=seurat_levels)

TRACK_ORDER <- TIS_SC_COMBO[order(df_NAMES$TISSUE, df_NAMES$SEURAT_CLUSTERS)]


cut_off <- dat_circcos$freq[length(dat_circcos$freq)*0.05]
TRACK_COLORS <- rep("black", nrow(dat_circcos))
TRACK_COLORS[dat_circcos$freq < cut_off] = "gray90"   #other NULL colour "#00000000"
TRACK_COLORS[c(21, 30)] = "red"  #C: 5, 8, 13 #B: 21, 30,  #A:7, 25, 47

#CIRCOS
pdf(file = paste0(out_dir, "/donor", donor, "_circos_disseminateduniquelinks_renamedsubset_red.pdf"), 
    width = 7, height = 7) 

circos.clear()
par(mar = rep(.9,4), cex=0.9)
circos.par(start.degree = 80, gap.degree = 1)

chordDiagramFromDataFrame(dat_circcos,
                          grid.col=GRID.COL,
                          order=TRACK_ORDER,
                          col=TRACK_COLORS,
                          directional=0,
                          annotationTrack = "grid",
                          annotationTrackHeight = 0.04,
                          preAllocateTracks = 1, 
                          link.decreasing = T)

for(si in get.all.sector.index()) {
  siLabel = gsub("\\|", " ", si)
  tissue <- substring(siLabel,1, 3)
  cluster <- substring(siLabel, 5)
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  if(tissue == "APP") {
    circos.text(mean(xlim), ylim[1], cluster, facing = "clockwise", adj = c(0, 0),
                niceFacing = TRUE, cex = 0.8, col = "blue", sector.index = si, 
                track.index = 1)
  } else if (tissue == "MLN") {
    circos.text(mean(xlim), ylim[1], cluster, facing = "clockwise", adj = c(0, 0),
                niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, 
                track.index = 1)
  } else if (tissue == "SPL") {
    circos.text(mean(xlim), ylim[1], cluster, facing = "clockwise", adj = c(0, 0),
                niceFacing = TRUE, cex = 0.8, col = "red", sector.index = si, 
                track.index = 1)
  }
}

dev.off()
