#Figure 5 (and S4)

#See Cytomapper for details about how to set up Images/ Masks/ and SCE object. Histocat (https://bodenmillergroup.github.io/histoCAT/) was used for composite images.

#Libraries
library(cytomapper)
library(ExperimentHub)
library(SingleCellExperiment)
library(S4Vectors)

library(scater)
library(ggpubr)
library(dplyr)
library(scran)
library(tidyr)

library(lme4)

#Directories
dir <- file.path()
out_dir <- file.path("figure images", "FigIMC")

##Analysis
spleenImages <- readRDS(file.path(dir, "spleen_images.rds"))
spleenMasks <- readRDS(file.path(dir, "spleen_masks.rds"))
spleenSCE <- readRDS(file.path(dir, "spleen_sce.rds"))

#Take out SPL_11
remove_mask <- match("SPL_11-Cell_Mask", names(spleenMasks))
spleenMasks <- getImages(spleenMasks, -remove_mask)

remove_image <- match("SPL_11-stack", names(spleenImages))
spleenImages <- getImages(spleenImages, -remove_image)

spleenSCE <- subset(spleenSCE, , ImageName!="SPL_11")

#Gating CD1c cells
cd1c_dat <- data.frame("CD1c" = retrieveCellInfo(spleenSCE, by = "CD1c", search = c("assays"), exprs_values = "logexprs")$value)
cd20_dat <- data.frame("CD20" = retrieveCellInfo(spleenSCE, by = "CD20", search = c("assays"), exprs_values = "logexprs")$value)

image_dat <- data.frame("ImageName" = retrieveCellInfo(spleenSCE, by = "ImageName", exprs_values = "logexprs")$value)
id_dat <- data.frame("id" = retrieveCellInfo(spleenSCE, by = "id", exprs_values = "logexprs")$value)

retrieve_dat <- data.frame(id_dat, image_dat, cd1c_dat, cd20_dat)

threshold <- 1.2

retrieve_dat$mz_label <- retrieve_dat$CD1c > threshold
colLabels(spleenSCE) <- retrieve_dat$mz_label

mz_SCE <- subset(spleenSCE, , label==TRUE)

#Fig 5c: Setting threshold scatter plot
contour_plot <- ggcells(spleenSCE, mapping = aes_string(y = "CD1c",x = "CD20"), exprs_values = "logexprs") +
  geom_point(size = 0.5, colour = "grey")+
  geom_density_2d(binwidth = 0.1)+
  geom_hline(yintercept = threshold, colour = "red")+
  #labs(subtitle = paste0("CD1c Treshold: ", threshold)) +
  theme_pubr()
file_name <- paste0("ContourPlot_CD1cThreshold_all.png")
ggsave(filename = file.path(out_dir, file_name), plot = contour_plot, width = 10, height = 10, units = c("cm"))

#Fig 5d + Supp 4a: Overlay of metadata on segmentation masks
cur_SPL1 <- subset(spleenSCE, , ImageName == "SPL_1")
mask_SPL1 <- getImages(spleenMasks, 1)

plotCells(mask = spleenMasks, object = spleenSCE,
          cell_id = "CellNumber", img_id = "ImageName", colour_by = "CD1c",
          outline_by = "label", display = "single",
          save_plot = list(filename = file.path(out_dir, "CD1cThreshold_CD1cExprs",
                                                "ImagesThreshold_CD1c.png")))

plotCells(mask = mask_SPL1, object = cur_SPL1,
          cell_id = "CellNumber", img_id = "ImageName", colour_by = "CD1c",
          outline_by = "label", display = "single", image_title = NULL,
          save_plot = list(filename = file.path(out_dir,
                                                "ImagesThreshold_CD1c_nolabel.png")))

#Fig 5e: Boxplot DDX21
DDX21_df <- makePerCellDF(spleenSCE, features = "DDX21", exprs_values = "logexprs")
ddx21_median <- DDX21_df %>% group_by(ImageName, label) %>% summarise(median = median(DDX21))

## Shapiro-Wilk normality test for the differences
shapiro.test(ddx21_median$median) # => p-value = 0.5598

rename <- data.frame(new_label = c("FALSE" = "CD1c-", "TRUE" = "CD1c+"))
rename$label <- rownames(rename)

ddx21_median <- merge(ddx21_median, rename, by = "label")

plot <- ggpaired(ddx21_median, x = "new_label", y = "median", color = "new_label", line.color = "gray", line.size = 0.8, id = "ImageName", point.size = 2) + 
  stat_compare_means(paired = TRUE, method = "t.test", label = "p.signif", label.x = 1.4, size = 7, label.y = 1.55) +
  labs(y = "DDX21 Expression Intensity") +
  scale_color_manual(values=c("#CE3234", "#BC80B8")) +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size=14))
file_name <- paste0("DDX21_boxplot.png")
ggsave(filename = file.path(out_dir, file_name), plot = plot, width = 7, height = 10, units = c("cm"))

#Fig 5g-i + Supp Fig 4b
exprs_df <- makePerCellDF(mz_SCE, features = rownames(mz_SCE), exprs_values = "logexprs")

order <- c("SPL_1", "SPL_2", "SPL_3", "SPL_4", "SPL_5", "SPL_6", "SPL_7", "SPL_8",
           "SPL_9", "SPL_10")

exprs_df$ImageName <- factor(exprs_df$ImageName,levels = order)

##Determining centre
gc_centre <- exprs_df %>% 
  group_by(ImageName) %>% 
  summarise(centre_X = mean(Pos_X),
            centre_Y = mean(Pos_Y))

exprs_geo_dat <- merge(exprs_df, gc_centre, by = "ImageName")

##Determining radius  
exprs_geo_dat <- exprs_geo_dat %>% group_by(ImageName) %>% mutate(r = ((Pos_X - centre_X)^2+(Pos_Y - centre_Y)^2)^0.5)

#Take out 5 and 7
exprs_geo_dat <- subset(exprs_geo_dat, ImageName != "SPL_5" & ImageName != "SPL_7")

#Scale radius by ImageName (max value = 1)
exprs_geo_dat <- exprs_geo_dat %>% group_by(ImageName) %>% mutate(norm_r = (r-min(r))/(max(r)-min(r)))

#scale such that 0.5 is mean and 1 = radius of image
exprs_geo_dat <- exprs_geo_dat %>% group_by(ImageName) %>% mutate(mean_r = mean(norm_r), norm_r2 = (0.5/mean_r)*norm_r)

ggplot(exprs_geo_dat, aes(x = Pos_X, y = Pos_Y)) +
  geom_point(size = 1) +
  geom_point(aes(x=centre_X, y=centre_Y), colour="red", size = 4) +
  facet_wrap(~ImageName) +
  theme_pubr() +
  labs(x = "X Position", y = "Y Position") +
  theme(axis.text = element_blank()) +
  scale_y_reverse()

library(rlang)
slope_dat <- {} 

for (i in 7:8){
  marker <- colnames(exprs_geo_dat)[i]
  
  test <- exprs_geo_dat %>% group_by(ImageName) %>% 
    mutate(temp_dat = !!sym(marker),
           fill = cut(temp_dat, 2, labels = c("Low", "High")))
  
  radial_plot <- ggplot(test, aes(x = norm_r2)) +
    geom_histogram(binwidth = 0.05, aes(fill = fill), position="dodge", show.legend = TRUE) +
    facet_wrap(~ImageName) +
    theme_pubr() +
    aes(y = ..ncount..) +
    labs(y = "Count (scaled)", x = "Radius (normalised)", fill = paste0(marker))
  
  file_name <- paste0("MZ_", marker, "proportion.png")
  ggsave(filename = file.path(out_dir, "Radial_MarkerExpression", file_name), plot = radial_plot, width = 30, height = 20, units = c("cm"))
  
  #
  cols_gg <- c("SPL_1" = "#00B0F6", "SPL_2" = "#00BF7D", "SPL_3" = "#00BFC4", "SPL_4" = "#39B600", "SPL_5" = "#9590FF", "SPL_6" = "#A3A500", "SPL_7" = "#D89000", "SPL_8" = "#E76BF3", "SPL_9" = "#F8766D", "SPL_10" = "#FF62BC") 
  
  radial_plot <- ggplot(test, aes(x = norm_r2, colour = ImageName)) +
    geom_histogram(binwidth = 0.06, aes(fill = fill), position="dodge", show.legend = TRUE) +
    scale_color_manual(values=cols_gg) +
    facet_wrap(~ImageName) +
    theme_pubr() +
    labs(y = "Count (scaled)", x = "Radius", fill = paste0(marker)) +
    aes(y = ..ncount..)
  
  pg <- layer_data(radial_plot) #x, xmin, xmax determines what bin cell is in
  
  stack_dat <- data.frame("cols" = pg$colour, "fill" = pg$fill, x = pg$x, y = pg$y, xmax = pg$xmax, xmin = pg$xmin, count = pg$ncount)
  
  cols <- data.frame("cols" = cols_gg)
  cols$ImageName <- rownames(cols)
  stack_dat <- merge(stack_dat, cols, by = c("cols"))
  
  exprs_string <- data.frame("Expression" = c("High", "Low"), "fill" = c("#00BFC4", "#F8766D"))
  stack_dat <- merge(stack_dat, exprs_string, by = c("fill"))
  stack_dat <- stack_dat[order(stack_dat$x),]
  
  prop_dat <- stack_dat %>% group_by(ImageName, Expression) %>% mutate("Bin" = rank(x))
  prop_dat <- prop_dat %>% group_by(ImageName, Bin) %>% mutate("Proportion" = y/sum(y))
  
  order <- c("SPL_1", "SPL_2", "SPL_3", "SPL_4", "SPL_5", "SPL_6", "SPL_7", "SPL_8",
             "SPL_9", "SPL_10")
  
  prop_dat$ImageName <- factor(prop_dat$ImageName,levels = order)
  
  prop_plot <- ggplot(prop_dat, aes(fill=Expression, y=y, x=Bin, width=1)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values = c("High" = "#00BFC4", "Low" = "#F8766D")) +
    facet_wrap(~ImageName) +
    theme_pubr() +
    labs(x = "Radius (normalised)", y = "Proportion", fill = paste0(marker)) +
    theme(axis.text.x=element_blank())
  
  file_name <- paste0("MZ_", marker, "proportion_total.png")
  ggsave(filename = file.path(out_dir, "Radial_MarkerExpression", file_name), plot = prop_plot, width = 30, height = 20, units = c("cm"))
  
  
  #Comparing slopes
  lm_dat <- subset(prop_dat, Expression == "High" & Proportion != 1)
  
  fits <- lmList(Proportion ~ x | ImageName, weights = count, data=lm_dat)
  
  slope_plot <- ggplot(lm_dat, aes(y=Proportion, x=x)) + 
    geom_point() +
    facet_wrap(~ImageName) +
    theme_pubr() +
    geom_smooth(method = "lm", mapping = aes(weight = count)) +
    labs(x = "Radius (normalised)", y = "Proportion (of high expression)", subtitle = marker, caption = "Linear regression, weighted by counts")
  
  file_name <- paste0("MZ_", marker, "proportion_slope.png")
  ggsave(filename = file.path(out_dir, "Radial_MarkerExpression", file_name), plot = slope_plot, width = 30, height = 20, units = c("cm"))
  
  #
  slope_length <- {}
  
  for(j in 1:length(unique(prop_dat$ImageName))){
    image_j <- names(fits)[j]
    slope_j <- fits[[j]]$coefficients[2]
    
    slope_length <- c(slope_length, slope_j)
    
  }
  temp_dat <- data.frame("Slope" = slope_length, "ImageName" = names(fits), "Marker" = marker)
  
  slope_dat <- rbind(slope_dat, temp_dat)
  
}

#Fig 5j: One sample t test DDX21 & DNA
marker_use_list <- unique(slope_dat$Marker)

p_values <- {}
marker_list <- {}

for(i in 1:length(unique(slope_dat$Marker))){
  marker <- marker_use_list[i]
  temp <- slope_dat[slope_dat$Marker == marker,]
  print(shapiro.test(temp$Slope)) #normal distribution
  
  res <- t.test(temp$Slope, mu = 0, alternative = "two.sided", conf.level = 0.95)
  
  p_values <- c(p_values, res$p.value)
  marker_list <- c(marker_list, marker)
}

ggplot(subset(slope_dat, Marker %in% c("DDX21", "DNA")), aes(x = Marker, y = Slope, fill = Marker)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", size = 1) +
  theme_pubr() +
  scale_fill_manual(values = c("#00BFC4", "#999999")) +
  labs(y = "Slope") +
  theme(axis.title.x = element_blank(), legend.position = "none")

#Fig 5K

##Dataframe with CD1c, CD20, DDX21 values
cd1c_dat <- data.frame("CD1c" = retrieveCellInfo(spleenSCE, by = "CD1c", search = c("assays"), exprs_values = "logexprs")$value)
cd20_dat <- data.frame("CD20" = retrieveCellInfo(spleenSCE, by = "CD20", search = c("assays"), exprs_values = "logexprs")$value)
ddx21_dat <- data.frame("DDX21" = retrieveCellInfo(spleenSCE, by = "DDX21", search = c("assays"), exprs_values = "logexprs")$value)

image_dat <- data.frame("ImageName" = retrieveCellInfo(spleenSCE, by = "ImageName", exprs_values = "logexprs")$value)
id_dat <- data.frame("id" = retrieveCellInfo(spleenSCE, by = "id", exprs_values = "logexprs")$value)

retrieve_dat <- data.frame(id_dat, image_dat, cd1c_dat, cd20_dat, ddx21_dat)

cd1c_threshold <- 1.2
cd20_threshold <- 1

retrieve_dat$bcell_label <- retrieve_dat$CD20 > cd20_threshold
retrieve_dat$mz_label <- retrieve_dat$CD1c > cd1c_threshold

mz_retrieve_dat <- subset(retrieve_dat, bcell_label == TRUE & mz_label == TRUE)

cut_quartile <- function(x){factor(findInterval(x, c(-Inf, quantile(x, probs=c(0.25, .5, .75)), Inf)), labels=c("DDX21_low", "DDX21_mid", "DDX21_mid","DDX21_high"))}

#Determine top / bottom X% quantile
mz_retrieve_dat <- mz_retrieve_dat %>% group_by(ImageName) %>% mutate(DDX21_quartile = cut_quartile(DDX21))

test <- left_join(retrieve_dat, mz_retrieve_dat)

label_mz <- test$DDX21_quartile

label_mz <- gtools::na.replace(as.vector(label_mz), "other")

colLabels(spleenSCE) <- label_mz

#Colour cells by label
plotCells(mask = spleenMasks, object = spleenSCE,
          cell_id = "CellNumber", img_id = "ImageName", 
          colour_by = "label", display = "single", 
          colour = list(label = c(DDX21_low = "blue", DDX21_high = "red", DDX21_mid = "grey", other = "grey")), image_title = NULL, scale_bar = list(length = 140),
          save_plot = list(filename = file.path(proj_dir, out_dir, "MZ_DDX21quartilesoverlay_imagetitlenull.png")))
