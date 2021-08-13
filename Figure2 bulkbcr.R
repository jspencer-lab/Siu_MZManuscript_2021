#Figure 2 bulk BCR

#Purpose: combine tissue/ datatype files, and graph results

#Resources
# - 5'RACE https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5879793/ 
# - Genewiz: https://www.genewiz.com/en-GB/Public/Services/Next-Generation-Sequencing/Immunogenomics/Immuno-Profiling
# - FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# - Trimmomatic manual: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

#Library
library(plyr) # plyr is needed for circo plot but need to be load before dplyr
library(tidyverse)
library(alakazam)
library(shazam)
library(cowplot)
library(rstatix)
library(ggpubr)
library(kableExtra)
library(glue)
library(airr)

library(RColorBrewer)
library(data.table)
library(foreach)

#Directories 
dat_dir <- file.path("Genewiz output", "Immcantation pipeline", "data", "results")
out_dir <- file.path("output")

#Combine tissues
donor <- "B" 

db_APP <- read_airr(file.path(dat_dir, paste0("22072021_", donor, "_Appendix_unique_db-pass_parse-select.tsv")))
db_SPL <- read_airr(file.path(dat_dir, paste0("22072021_", donor, "_Spleen_unique_db-pass_parse-select.tsv")))

db_all <- dplyr::bind_rows(list(APP = db_APP, SPL = db_SPL), .id = 'tissue')

db_all <- db_all %>% dplyr::mutate(donor = paste0(donor2))

alakazam::writeChangeoDb(db_all, file = file.path(dat_dir, paste0("donor", donor, "-alltissues_parse-select.tsv")))

#Determine clone
donor <- "C"
dat <- readChangeoDb(file = file.path(dat_dir, paste0("donor", donor, "-alltissues_parse-select.tsv")))
db <- distToNearest(dat, model="ham", normalize="len", vCallColumn="v_call", nproc=16)

# Determine threshold
threshold <- findThreshold(db$dist_nearest, method="density")
thr <- round(threshold@threshold, 2)

p1 <- ggplot(subset(db, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.15, color="blue", linetype=2) +
  #geom_vline(xintercept=0.10, color="cyan", linetype=2) +
  labs(subtitle = paste0("Donor ", donor), caption = paste0("Calculated threshold = ", thr))
plot(p1)

ggsave(filename = file.path(out_dir, paste0("ThresholdPlot_donor", donor, ".png")), plot = p1)

#Combine bulk and singlecell seq files for clonal definition
dat_dir <- file.path("Genewiz output", "Immcantation pipeline", "data", "results")
scVDJ_dir <- file.path("Projects", "2020_tissues10X")
seurat_metadatadir <- file.path("Projects", "2020_tissues10X", "Data", "Seurat saved objects")

comboout_dir <- file.path("Genewiz output", "Immcantation pipeline", "data", "results", "combo")

donor <- "C"
dat <- readChangeoDb(file = file.path(dat_dir, paste0("donor", donor, "_combined_clone-pass.tsv")))

db <- readChangeoDb(file.path(scVDJ_dir, paste0("donor", donor, "_data_ph_genotyped_germ-pass.tsv")))

#Rename cregion to c_call in bulk
dat$c_call <- sub("Human-", "", dat$cregion)
dat$c_call <- sub("-InternalC", "", dat$c_call)

#Rename bulk v/j_call to match 10x
temp <- str_match_all(dat$v_call, regex("IGHV[0-9]-[0-9]*\\*[0-9]*"))
temp <- t(as.data.frame(lapply(temp, FUN=function(w){
  paste(w[,1], collapse = ', ')
})))
rownames(temp) <- rownames(dat)
dat <- subset(dat, select = -c(v_call))
colnames(temp) <- 'v_call'
dat <- cbind(dat, temp)

temp2 <- str_match_all(dat$j_call, regex("IGHJ[0-9]*\\*[0-9]*"))
temp2 <- t(as.data.frame(lapply(temp2, FUN=function(w){
  paste(w[,1], collapse = ', ')
})))
rownames(temp2) <- rownames(dat)
dat <- subset(dat, select = -c(j_call))
colnames(temp2) <- 'j_call'
dat <- cbind(dat, temp2)

#Rename clone id to bulk_cloneid vs sc_cloneid
names(dat)[names(dat) == 'clone_id'] <- 'bulk_cloneid'
names(db)[names(db) == 'clone_id'] <- 'sc_cloneid'

#Combine
combo_dat <- dplyr::bind_rows(list(bulk = dat, sc = db), .id = 'data_type')

alakazam::writeChangeoDb(combo_dat, file = file.path(comboout_dir, paste0("donor", donor, "-combo_dat.tsv")))

##COMBO DISTANCE
combo_dat_distance <- distToNearest(combo_dat, model="ham", normalize="len", vCallColumn="v_call", nproc=16)

# Determine threshold
threshold <- findThreshold(combo_dat_distance$dist_nearest, method="density")
thr <- round(threshold@threshold, 2)

p1 <- ggplot(subset(combo_dat_distance, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.15, color="blue", linetype=2) +
  #geom_vline(xintercept=0.10, color="cyan", linetype=2) +
  labs(subtitle = paste0("Donor ", donor), caption = paste0("Calculated threshold = ", thr))
plot(p1)

ggsave(filename = file.path(comboout_dir, paste0("ComboThresholdPlot_donor", donor, ".png")), plot = p1)

#Fig 2b and S3G
df_all <- {}
donor_list <- c("A", "B", "C")
version <- "version1"

for (i in 1:length(donor_list)){
  donor <- donor_list[i]
  dat <- readChangeoDb(file = file.path(dat_dir, version, paste0("donor", donor, "_datcombined_clone-pass.tsv")))
  
  #Rename IGHA1/2 to just IGHA etc
  dat$isotypes <- forcats::fct_collapse(dat$c_call, IGHG = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), IGHA = c("IGHA1", "IGHA2"), IGLC = c("IGLC-1", "IGLC-2", "IGLC-3", "IGLC-4", "IGLC1"))
  
  #Adding cell phenotypes into data where possible via cell_id
  seurat_metadata <- readRDS(file.path(seurat_metadatadir, paste0("integratedall_metadata_donor", donor, ".rds")))
  
  seurat_metadata  <- rownames_to_column(seurat_metadata , "cell_id")
  seurat_metadata <- seurat_metadata[!duplicated(seurat_metadata$cell_id), ]
  
  subset_ids <- data.frame(seurat_clusters = 0:22, 
                           subset_names = c("MZB-1", "TS", "ABC3", "ABC1", "Memory", "aNAV", "Naive", "Memory", "ABC2", "IgM-only", "MZB-2", "Naive", "Memory", "GC", "IgM-only", "GC", "Naive", "GC", "DN-A", "aNAV", "DN-B", "ABC4", "PB"))
  
  subset_ids$seurat_clusters <- as.factor(subset_ids$seurat_clusters)
  
  seurat_metadata <- dplyr::right_join(seurat_metadata, subset_ids, by = c("seurat_clusters"))
  
  add_dat <- data.frame("subset_names" = seurat_metadata$subset_names, "cell_id" = seurat_metadata$cell_id, "seurat_clusters" = seurat_metadata$seurat_clusters)
  
  dat <- left_join(dat, add_dat, by='cell_id')
  
  dat <- subset(dat, tissue != "MLN")
  
  db <- dat[!(dat$subset_names == "MZB-1" & dat$isotypes == "IGHA") & 
              !(dat$subset_names == "MZB-1" & dat$isotypes == "IGHG") &
              !(dat$subset_names == "MZB-2" & dat$isotypes == "IGHA") &
              !(dat$subset_names == "MZB-2" & dat$isotypes == "IGHG"),]
  
  dt <- data.table(db)
  groupByFormatted <- 'clone_id'
  dt <- dt[, list(yidx = list(.I)) , by = groupByFormatted] #yidx = indices of repeats
  groups <- dt[, yidx]
  groupsLength <- unlist(lapply(groups,length))
  singletons <- groups[groupsLength==1]
  groups <- groups[groupsLength>1]
  numbOfTotalGroups <- length(groups)
  
  listResults  <- foreach(i = iterators::icount(numbOfTotalGroups)) %do% {
    dfClone <- db[groups[[i]], c("cdr3", "clone_id", "sc_cloneid", "bulk_cloneid", "tissue", "isotypes", "data_type", "subset_names", "seurat_clusters")] 
    
    dfClone <- unique(dfClone)
    
    matReturn <- NULL
    temp_datatype <- NULL
    temp_scclone <- NULL
    temp_bulkclone <- NULL
    temp_clone <- NULL
    temp_res <- NULL
    datatype <- NULL
    cloneid_10x <- NULL
    
    #Is there both datasets?
    temp_datatype <- unique(dfClone$data_type)
    temp_scclone <- unique(dfClone$sc_cloneid)
    temp_bulkclone <- unique(dfClone$bulk_cloneid)
    
    if (length(temp_datatype) > 1){
      datatype <- "both"
      temp_scclone <- temp_scclone[!is.na(temp_scclone)]
      temp_bulkclone <- temp_bulkclone[!is.na(temp_bulkclone)]
    } else if (length(temp_datatype) ==1){
      datatype <- temp_datatype
    }
    
    #Is there more than one sc clone
    if (length(temp_scclone) > 1){
      cloneid_10x <- "morethan1"
    } else if (length(temp_scclone) ==1){
      cloneid_10x <- temp_scclone
    }
    #Is there one bulk_cloneid?
    if (length(temp_bulkclone) > 1){
      cloneid_bulk <- "morethan1"
    } else if (length(temp_bulkclone) ==1){
      cloneid_bulk <- temp_bulkclone
    }
    
    #Is it resident or disseminate? And what type of dissemination
    temp_res <- unique(dfClone$tissue)
    if (length(temp_res) > 1){
      clonetype <- "disseminate"
      
      tissue_bulk <- unique(subset(dfClone, data_type == "bulk", select = "tissue"))
      tissue_sc <- unique(subset(dfClone, data_type == "sc", select = "tissue"))
      
      #Determining score
      score <- 0
      if (nrow(tissue_bulk) > 1){
        score <- score + 1
      }
      if (nrow(tissue_sc) > 1){
        score <- score + 10
      }
      if(nrow(tissue_sc) != 0 & nrow(tissue_bulk) != 0){
        if (!(tissue_bulk %in% tissue_sc)){
          score <- score + 100
        }
      }
    }else if (length(temp_res) ==1){
      clonetype <- "resident"
      score <- "NA"
    }
    
    
    #Linking to MZB/ subsets
    if ("MZB-1" %in% dfClone$subset_names & "MZB-2" %in% dfClone$subset_names){
      MZB_link <- "both"
    } else if ("MZB-1" %in% dfClone$subset_names){
      MZB_link <- "MZB-1"
    } else if("MZB-2" %in% dfClone$subset_names){
      MZB_link <- "MZB-2"
    } else {
      MZB_link <- "NA"
    }
    
    #String of associated subsets
    subset_string <- sort(unique(dfClone$subset_names[!is.na(dfClone$subset_names)]))
    subset_string <- paste(subset_string, collapse=",")
    
    matReturn <- cbind("clone_id" = unique(dfClone$clone_id), cloneid_10x, cloneid_bulk, clonetype, score, datatype, MZB_link, donor, subset_string)
    return(matReturn)
  }
  
  df <- do.call(rbind, listResults)
  
  #COMBINE ALL DONORS
  df_all <- rbind(df_all, df)
}

df_all <- data.frame(df_all)

#Converting label
score_convert <- data.frame(score = c(NA, 1, 10, 11, 100, 101, 110, 111),
                            disseminationtype = c("none", "bulk", "sc", "bulk_sc", "diff", "diff_bulk", "diff_sc", "diff_bulk_sc"))
score_convert <- score_convert %>%  
  mutate(score = as.character(score))

df_all <- dplyr::left_join(df_all, score_convert, by = "score")


#Final stats with separate clonal definitions
df_all2 <- subset(df_all, MZB_link %in% c("MZB-1", "MZB-2"))

#Categories
df_all2$bulk <- forcats::fct_collapse(df_all2$disseminationtype, bulk = c("bulk_sc", "diff_bulk", "diff_bulk_sc"), resident = c("diff", "diff_sc", "sc"))

df_all2$sc <- forcats::fct_collapse(df_all2$disseminationtype, sc = c("bulk_sc", "diff_sc", "diff_bulk_sc"), resident = c("diff", "diff_bulk", "bulk"))

df_all2$diff <- forcats::fct_collapse(df_all2$disseminationtype, diff = c("diff", "diff_bulk", "diff_sc", "diff_bulk_sc"), resident = c("bulk_sc", "sc", "bulk"))

#replace NA with other in categories
df_all2$bulk[is.na(df_all2$bulk)] <- "resident"
df_all2$sc[is.na(df_all2$sc)] <- "resident"
df_all2$diff[is.na(df_all2$diff)] <- "resident"

df_all2$disseminate <- df_all2$clonetype
df_all3 <- subset(df_all2, datatype == "both")

#Loop
all_savestat <- {}
dis_types <- c("bulk", "sc", "diff", "disseminate")

for (i in 1:length(dis_types)){
  type <- dis_types[i]
  
  sum_both <- df_all3 %>% group_by(donor, MZB_link, !!sym(type)) %>% count()
  sum_all <- df_all2 %>% group_by(donor, MZB_link, !!sym(type)) %>% count()
  
  if (nrow(sum_both) < 12){
    sum_both <- sum_both %>% ungroup() %>% add_row(donor = "B", MZB_link = "MZB-2", sc = "sc", n = 0) %>% group_by(donor, MZB_link, !!sym(type))
  }
  
  if (nrow(sum_all) < 12){
    sum_all <- sum_all %>% ungroup() %>% add_row(donor = "B", MZB_link = "MZB-2", sc = "sc", n = 0) %>% group_by(donor, MZB_link, !!sym(type))
  }
  
  both_stat <- sum_both %>% group_by(donor, MZB_link) %>% summarise(dis = n[!!sym(type) == paste0(type)], res = n[!!sym(type) == "resident"]) %>% mutate(ratio = dis/ res)
  
  all_stat <- sum_all %>% group_by(donor, MZB_link) %>% summarise(dis = n[!!sym(type) == paste0(type)], res = n[!!sym(type) == "resident"]) %>% mutate(ratio = dis/ res)
  
  both_stat$id <- "both"
  all_stat$id <- "all"
  
  save_stat <- rbind(both_stat, all_stat)
  save_stat$type <- paste0(type)
  
  all_savestat <- rbind(all_savestat, save_stat)
}

test <- all_savestat %>% group_by(donor, id, type) %>% summarise(mzb1_dis = dis[MZB_link == "MZB-1"], mzb1_res = res[MZB_link == "MZB-1"], mzb2_dis = dis[MZB_link == "MZB-2"], mzb2_res = res[MZB_link == "MZB-2"], mzb1_ratio = ratio[MZB_link == "MZB-1"], mzb2_ratio = ratio[MZB_link == "MZB-2"]) %>% mutate(ratiomzb2_mzb1 = mzb2_ratio/ mzb1_ratio)

write.csv(test, file = file.path(out_dir, version, paste0(Sys.Date(),"_dissemination_ratios_onlyappspl.csv")))

plot_dat <- subset(test, id == "all")

plot <- plot_dat %>% mutate(type = factor(type, levels = c("sc", "diff", "bulk", "disseminate"))) %>%
  ggplot(aes(x = type, y = ratiomzb2_mzb1, col = donor)) +
  geom_point(size = 5) +
  theme_pubr() +
  ylim(-0.01, 1.2) +
  labs(x = "Estimations of relative dissemination of MZB1 and MZB2 \n defined in SC data", y = "Dissemination ratio \n MZB-2 / MZB-1", col = "Donor") +
  theme(axis.title.x = element_text(face="bold")) +
  geom_hline(yintercept = 1, linetype='dotted') +
  scale_x_discrete(labels = c('Clones in \n both tissues \n of SC only', 'SC with \n relatives in the \n other tissue \n in bulk data','SC with \n relatives in \n both tissues \n in bulk data', 'Combined \n dissemination \n index'))

ggsave(filename = file.path(out_dir, version, paste0(version, "plot_disseminationratio.png")), plot = plot, width = 15, height = 10, units = "cm")
