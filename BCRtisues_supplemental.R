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
out_dir <- file.path("C:", "Users", "jacqu", "Dropbox", "NI Submission 2021", "figure images", "Supplemental", "FigS5")

##Set up
db_all_isotype <- subset(db_all, c_call %in% c("IGHM","IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2"))

#Rename isotypes
db_all_isotype$isotype <- forcats::fct_collapse(db_all_isotype$c_call, IGHG = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), IGHA = c("IGHA1", "IGHA2"))


##Abundance
groupby <- "isotype"
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
groupby <- "subset_names"
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