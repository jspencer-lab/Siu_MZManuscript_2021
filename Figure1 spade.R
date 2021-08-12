#Figure1 spade

#Libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

#Directories
out_dir <- file.path("figure images", "Revisions")
in_dir <- file.path("Projects", "tissuesCyTOF")

#Import (data exported as csv from cytobank)
dat <- read.csv(file = file.path(in_dir, "Cytobank_exportedbubbles_exportedstats.csv"))

#Figure 1H
cd19_counts <- dat %>% subset(Conditions == "CD19+") %>% group_by(Donor, Tissues) %>% select(Cell_counts)
colnames(cd19_counts) <- c("Donor", "Tissues", "CD19_counts")

dat <- merge(dat, cd19_counts, by = c("Donor", "Tissues")) 

dat <- dat %>% group_by(Donor, Tissues) %>% mutate(CD19_prop = Cell_counts / CD19_counts * 100)

plot_dat <- subset(dat, Conditions == "MZB1"| Conditions == "MZB2")

bxp <- ggboxplot(plot_dat, x = "Conditions", y = "CD19_prop",
                 color = "Tissues", palette = c("Appendix" = "blue","mLN" = "black","Spleen" = "red"), add = "jitter", shape = "Tissues", outlier.shape = NA
) + labs(x = "Subset Names", y = "% of CD19", color = "Tissues", shape = "Tissues") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x=element_blank(), legend.position = "right")
file_name <- paste0("CyTOFSPADE_MZBAbundance.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = bxp, width = 12, height = 12, units = c("cm"))

#Figure 1G
mz_dat <- subset(dat, Conditions == "MZB1" | Conditions == "MZB2")
mz_dat$sampleid <- paste0(mz_dat$Tissues, "_", mz_dat$Donor)

markers <- colnames(mz_dat)[c(8, 10:11)]

p <- list()
for (i in 1:length(markers)){
  marker <- markers[[i]]
  p[[i]] <- ggpaired(mz_dat, x = "Conditions", y = marker, id = "sampleid",
                     color = "Conditions", line.color = "gray", line.size = 0.4, palette = c("MZB1" = "#009681","MZB2" = "#6CDC52"), shape = 21)+
    stat_compare_means(paired = TRUE, hide.ns = FALSE, method = "wilcox.test", label.x.npc = "centre", aes(label = ..p.signif..)) +
    labs(y = paste0("Median ", marker, " Expression")) +
    theme(legend.position = "none", axis.title.x=element_blank()) + 
    geom_point(aes_string(shape = "Conditions", fill = "Conditions", size = 3)) +
    scale_shape_manual(values=c(22, 24))
}

save_plot <- cowplot::plot_grid(plotlist = p, nrow = 1)

file_name <- paste0("MZ_MedianMarkerExpression_Cytobank.pdf")
ggsave(filename = file.path(out_dir, file_name), plot = save_plot, width = 20, height = 10, units = c("cm"))