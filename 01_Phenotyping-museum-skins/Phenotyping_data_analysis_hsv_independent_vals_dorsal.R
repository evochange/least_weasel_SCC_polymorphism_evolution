## Phenotyping dataset - museum skins
## PCA/KW analyses - HSV data
## Last edited: 12-05-2024, Inês Miranda


## Set working directory:
setwd("/pth/to/folder")
rm(list=ls())


## Import libraries:
library(factoextra)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(rstatix)
library(tidyverse)
library(ggpubr)

#### 1. TOTAL DATASET ANALYSIS ####

## Read input:
m <- read.table("HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt", 
                header=T, sep = "\t", na.strings = "n/a")
str (m)

## Check columns for amount of missing data:
## No missing data left because rgb2hsv was applied to filtered dataset already
hsv_ind_total <- m[complete.cases(m), ]
str(hsv_ind_total)
summary (hsv_ind_total)

## Info on number of specimens per category:
table(hsv_ind_total$SEASON,hsv_ind_total$MORPH)

## Prepare matrix for PCA analysis:
hsv_ind_M <- as.matrix (hsv_ind_total[,c(5:16)])
str(hsv_ind_M)

## Run PCA analysis:
rgb_pca <- prcomp (hsv_ind_M, scale. = T)

## Check PCA components:
summary(rgb_pca)
rgb_pca$x

## Extract percentage of variance values from PCA summary:
perc_variance <- as.vector (summary(rgb_pca)$importance[2,])
str(perc_variance)

## Check contribution of each variable to each PC:

  ## Get PCA variables:
  rbg_vars <- get_pca_var(rgb_pca)
  rbg_vars

  ## Check variables contributions to PCs:
  fviz_contrib(rgb_pca, choice = "var", axes = 1, top = 30)
  fviz_contrib(rgb_pca, choice = "var", axes = 2, top = 30)

## Create dataframe with scores:
scores = as.data.frame(rgb_pca$x)
str(scores)

indv = hsv_ind_total[,c(1:4)]
str (indv)

scores_indv = cbind (indv,scores)
str(scores_indv)


## Plot results:
ggplot(scores_indv, aes(x=PC1, y=PC2, shape=Group, fill=Group)) +
  theme_bw() +
  geom_point(size = 6) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_manual(values=c("#AD7135","#ffffff","#593A1B", "#cdb79e")) +
  xlab(paste0("PC1 (", round(perc_variance[1]*100,2),"%)")) +
  ylab(paste0("PC2 (", round(perc_variance[2]*100,2),"%)"))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 14)) +
  guides (shape = guide_legend (order=1), fill = guide_legend(order=2)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


#### 2. BROWN MORPHS DATASET ####
#rm(list=ls())

#Read input: [if starting code here]
#m <- read.table("HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt", 
#                header=T, sep = "\t", na.strings = "n/a")
#str (m)

## Subset data matrix:
m_brown <- m[m$Group != "nivalis winter",]
str(m_brown)

## Check columns for amount of missing data:
## No missing data left because rgb2hsv was applied to filtered dataset already
hsv_brown <- m_brown[complete.cases(m_brown), ]

str(hsv_brown)
summary (hsv_brown)

## Info on number of specimens per category:
table(hsv_brown$SEASON,hsv_brown$MORPH)


##### PLOT DATA #####

## Prepare matrix for PCA analysis:
hsv_brown_M <- as.matrix (hsv_brown[,c(5:16)])
str(hsv_brown_M)

## Run PCA analysis:
rgb_brown_pca <- prcomp (hsv_brown_M, scale. = T)

## Check PCA components:
summary(rgb_brown_pca)
rgb_brown_pca$x

## Correlation between original variables and PC components:
cor(hsv_brown_M, rgb_brown_pca$x)

## Extract percentage of variance values from PCA summary:
perc_brown_variance <- as.vector (summary(rgb_brown_pca)$importance[2,])
str(perc_brown_variance)

## Check contribution of each variable to each PC:

  ## Get PCA variables:
  rbg_brown_vars <- get_pca_var(rgb_brown_pca)
  rbg_brown_vars

  ## Check variables contributions to PCs:
  fviz_contrib(rgb_brown_pca, choice = "var", axes = 1, top = 30)
  fviz_contrib(rgb_brown_pca, choice = "var", axes = 2, top = 30)

## Create dataframe with scores:
scores_brown = as.data.frame(rgb_brown_pca$x)
str(scores_brown)

indv_brown = hsv_brown[,c(1:4)]
str (indv_brown)

scores_indv_brown = cbind (indv_brown,scores_brown)
str(scores_indv_brown)

## Plot results:
ggplot(scores_indv_brown, aes(x=PC1, y=PC2, shape=Group, fill=Group)) +
  theme_bw() +
  geom_point(size = 6) +
  scale_shape_manual(values = c(21, 23, 24)) +
  scale_fill_manual(values=c("#AD7135","#593A1B", "#cdb79e")) +
  xlab(paste0("PC1 (", round(perc_brown_variance[1]*100,2),"%)")) +
  ylab(paste0("PC2 (", round(perc_brown_variance[2]*100,2),"%)"))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 14)) +
  guides (shape = guide_legend (order=1), fill = guide_legend(order=2)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


## Biplot of multivariate data:

  ## Add group columns:
  brown_groups <- as.factor (paste (scores_indv_brown$MORPH,scores_indv_brown$SEASON, sep = " "))

  ## Plot Biplot:
  ggbiplot(rgb_brown_pca, choices = 1:2, groups = brown_groups, varname.size = 4, labels.size = 4, ellipse = TRUE,
          ellipse.alpha = 0) +
    theme_bw() +
    scale_color_manual(name = "Phenotypic group", values=c("#AD7135","#593A1B", "#cdb79e")) +
    scale_fill_manual(name = "Phenotypic group", values=c("#AD7135","#593A1B", "#cdb79e")) +
    scale_shape_manual(name = "Phenotypic group", values = c(21, 23, 24)) +
    geom_point(aes(color= brown_groups, fill=brown_groups, shape=brown_groups), size = 4) +
    xlab(paste0("PC1 (", round(perc_brown_variance[1]*100,2),"%)")) +
    ylab(paste0("PC2 (", round(perc_brown_variance[2]*100,2),"%)")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size = 12)) +
    xlim (-2.5, 2.5)
  


##### TEST ASSUMPTIONS FOR PAIRWISE ANOVA TESTS #####

## Check for univariate outliers:
hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(H_h)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(H_s)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(H_v)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(UD_h)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(UD_s)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(UD_v)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(MD_h)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(MD_s)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(MD_v)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(LD_h)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(LD_s)

hsv_brown %>%
  group_by(MORPH, SEASON) %>%
  identify_outliers(LD_v)


## Check univariate normality assumption:
hsv_brown %>%
  group_by(Group) %>%
  shapiro_test(H_h, H_s, H_v) %>%
  arrange(variable)
# VARIABLE H_h has a small deviation

hsv_brown %>%
  group_by(Group) %>%
  shapiro_test(UD_h, UD_s, UD_v) %>%
  arrange(variable)

hsv_brown %>%
  group_by(Group) %>%
  shapiro_test(MD_h, MD_s, MD_v) %>%
  arrange(variable)
# VARIABLE MD_S has a small deviation

hsv_brown %>%
  group_by(Group) %>%
  shapiro_test(LD_h, LD_s, LD_v) %>%
  arrange(variable)
# VARIABLE LD_S has a small deviation

## QQ plot of each variable with P<0.05:
ggqqplot(hsv_brown, "H_h", facet.by = "Group",
         ylab = "", ggtheme = theme_bw())

ggqqplot(hsv_brown, "MD_s", facet.by = "Group",
         ylab = "", ggtheme = theme_bw())

ggqqplot(hsv_brown, "LD_s", facet.by = "Group",
         ylab = "", ggtheme = theme_bw())



#### 3. VULGARIS MORPHS DATASET ####
#rm(list=ls())

#Read input:  [if starting code here]
#m <- read.table("HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt", 
#                header=T, sep = "\t", na.strings = "n/a")
#str (m)

## Subset data matrix:
m_vulg <- m[m$MORPH != "nivalis",]
str(m_vulg)

## Check columns for amount of missing data:
## No missing data left because rgb2hsv was applied to filtered dataset already
hsv_vulg <- m_vulg[complete.cases(m_vulg), ]

str(hsv_vulg)
summary (hsv_vulg)

## Info on number of specimens per category:
table(hsv_vulg$SEASON,hsv_vulg$MORPH)

## Prepare matrix for PCA analysis:
hsv_vulg_M <- as.matrix (hsv_vulg[,c(5:16)])
str(hsv_vulg_M)

## Run PCA analysis:
rgb_vulg_pca <- prcomp (hsv_vulg_M, scale. = T)

## Check PCA components:
summary(rgb_vulg_pca)
rgb_vulg_pca$x

## Extract percentage of variance values from PCA summary:
perc_vulg_variance <- as.vector (summary(rgb_vulg_pca)$importance[2,])
str(perc_vulg_variance)

## Check contribution of each variable to each PC:

  ## Get PCA variables:
  rbg_vulg_vars <- get_pca_var(rgb_vulg_pca)
  rbg_vulg_vars
  
  ## Check variables contributions to PCs:
  fviz_contrib(rgb_vulg_pca, choice = "var", axes = 1, top = 30)
  fviz_contrib(rgb_vulg_pca, choice = "var", axes = 2, top = 30)

## Create dataframe with scores:
scores_vulg = as.data.frame(rgb_vulg_pca$x)
str(scores_vulg)

indv_vulg = hsv_vulg[,c(1:4)]
str (indv_vulg)

scores_indv_vulg = cbind (indv_vulg,scores_vulg)
str(scores_indv_vulg)

## Plot results:
ggplot(scores_indv_vulg, aes(x=PC1, y=PC2, fill=SEASON, shape=SEASON)) +
  theme_bw() +
  geom_point(size = 6) +
  scale_shape_manual(name = "Season", values = c(23, 24)) +
  scale_fill_manual(name = "Season", values=c("#593A1B", "#cdb79e")) +
  xlab(paste0("PC1 (", round(perc_vulg_variance[1]*100,2),"%)")) +
  ylab(paste0("PC2 (", round(perc_vulg_variance[2]*100,2),"%)"))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 14)) +
  guides (shape = guide_legend (order=1), fill = guide_legend(order=2)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) 


## Biplot of multivariate data:

  ## Add group columns
  vulg_groups <- as.factor (paste (scores_indv_vulg$MORPH,scores_indv_vulg$SEASON, sep = " "))

  ## Plot Biplot:
  ggbiplot(rgb_vulg_pca, choices = 1:2, groups = vulg_groups, varname.size = 4, labels.size = 4, ellipse = TRUE,
          ellipse.alpha = 0) +
    theme_bw() +
    scale_color_manual(name = "Phenotypic group", values=c("#593A1B","#cdb79e")) +
    scale_fill_manual(name = "Phenotypic group", values=c("#593A1B","#cdb79e")) +
    scale_shape_manual(name = "Phenotypic group", values = c(23, 24)) +
    geom_point(aes(color= vulg_groups, fill=vulg_groups, shape=vulg_groups), size = 4) +
    xlab(paste0("PC1 (", round(perc_vulg_variance[1]*100,2),"%)")) +
    ylab(paste0("PC2 (", round(perc_vulg_variance[2]*100,2),"%)")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size = 12)) +
    xlim (-2.5, 2.2)


## RUN KW TEST:

  ## Group the data by variable
  grouped.vulg <- hsv_vulg %>%
    gather(key = "variable", value = "value", H_h, UD_h, MD_h, LD_h, H_s, UD_s, MD_s, LD_s, H_v, UD_v, MD_v, LD_v) %>%
    group_by(variable)

  ## Do Kruskal-Wallis test
  grouped.vulg %>% kruskal_test(value ~ SEASON)


#### 4. SUMMER MORPHS DATASET ####
#rm(list=ls())

#Read input: [if starting code here]
#m <- read.table("HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt", 
#                header=T, sep = "\t", na.strings = "n/a")
#str (m)

## Subset data matrix:
m_sum <- m[m$SEASON == "summer",]
str(m_sum)

## Check columns for amount of missing data:
## No missing data left because rgb2hsv was applied to filtered dataset already
hsv_sum <- m_sum[complete.cases(m_sum), ]

str(hsv_sum)
summary (hsv_sum)

## Info on number of specimens per category:
table(hsv_sum$SEASON,hsv_sum$MORPH)

## Prepare matrix for PCA analysis:
hsv_sum_M <- as.matrix (hsv_sum[,c(5:16)])
str(hsv_sum_M)

## Run PCA analysis:
rgb_sum_pca <- prcomp (hsv_sum_M, scale. = T)

## Check PCA components:
summary(rgb_sum_pca)
rgb_sum_pca$x

## Extract percentage of variance values from PCA summary:
perc_sum_variance <- as.vector (summary(rgb_sum_pca)$importance[2,])
str(perc_sum_variance)

## Check contribution of each variable to each PC:

  ## Get PCA variables:
  rbg_sum_vars <- get_pca_var(rgb_sum_pca)
  rbg_sum_vars

  ## Check variables contributions to PCs:
  fviz_contrib(rgb_sum_pca, choice = "var", axes = 1, top = 30)
  fviz_contrib(rgb_sum_pca, choice = "var", axes = 2, top = 30)

## Create dataframe with scores:
scores_sum = as.data.frame(rgb_sum_pca$x)
str(scores_sum)

indv_sum = hsv_sum[,c(1:4)]
str (indv_sum)

scores_indv_sum = cbind (indv_sum,scores_sum)
str(scores_indv_sum)

## Plot results:
ggplot(scores_indv_sum, aes(x=PC1, y=PC2, fill=MORPH, shape=MORPH)) +
  theme_bw() +
  geom_point(size = 6) +
  scale_shape_manual(name = "Season", values = c(21, 23)) +
  scale_fill_manual(name = "Morph", values=c("#AD7135","#593A1B")) +
  xlab(paste0("PC1 (", round(perc_sum_variance[1]*100,2),"%)")) +
  ylab(paste0("PC2 (", round(perc_sum_variance[2]*100,2),"%)"))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 14)) +
  guides (shape = guide_legend (order=1), fill = guide_legend(order=2)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


## Biplot of multivariate data:

  ## Add group columns:
  sum_groups <- as.factor (paste (scores_indv_sum$MORPH,scores_indv_sum$SEASON, sep = " "))

  ## Plot Biplot:
  ggbiplot(rgb_sum_pca, choices = 1:2, groups = sum_groups, varname.size = 4, labels.size = 4, ellipse = TRUE,
          ellipse.alpha = 0) +
    theme_bw() +
    scale_color_manual(name = "Phenotypic group", values=c("#AD7135","#593A1B")) +
    scale_fill_manual(name = "Phenotypic group", values=c("#AD7135","#593A1B")) +
    scale_shape_manual(name = "Phenotypic group", values = c(21, 23)) +
    geom_point(aes(color= sum_groups, fill=sum_groups, shape=sum_groups), size = 4) +
    xlab(paste0("PC1 (", round(perc_sum_variance[1]*100,2),"%)")) +
    ylab(paste0("PC2 (", round(perc_sum_variance[2]*100,2),"%)")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size = 12))
  
## RUN KW TEST:

  ## Group the data by variable
  grouped.sum <- hsv_sum %>%
    gather(key = "variable", value = "value", H_h, UD_h, MD_h, LD_h, H_s, UD_s, MD_s, LD_s, H_v, UD_v, MD_v, LD_v) %>%
    group_by(variable)

  ## Do Kruskal-Wallis test
  grouped.sum %>% kruskal_test(value ~ MORPH)



#### 5. OPPOSITE MORPHS DATASET ####
#rm(list=ls())

#Read input: [if starting code here]
#m <- read.table("HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt", 
#                header=T, sep = "\t", na.strings = "n/a")
#str (m)

## Subset data matrix:
m_oppo <- m[(m$Group == "nivalis summer") | (m$Group == "vulgaris winter"),]
str(m_oppo)

## Check columns for amount of missing data:
## No missing data left because rgb2hsv was applied to filtered dataset already
hsv_oppo <- m_oppo[complete.cases(m_oppo), ]

str(hsv_oppo)
summary (hsv_oppo)

## Info on number of specimens per category:
table(hsv_oppo$SEASON,hsv_oppo$MORPH)

## Prepare matrix for PCA analysis:
hsv_oppo_M <- as.matrix (hsv_oppo[,c(5:16)])
str(hsv_oppo_M)

## Run PCA analysis:
rgb_oppo_pca <- prcomp (hsv_oppo_M, scale. = T)

## Check PCA components:
summary(rgb_oppo_pca)
rgb_oppo_pca$x

## Extract percentage of variance values from PCA summary:
perc_oppo_variance <- as.vector (summary(rgb_oppo_pca)$importance[2,])
str(perc_oppo_variance)

## Check contribution of each variable to each PC:
  
  ## Get PCA variables:
  rbg_oppo_vars <- get_pca_var(rgb_oppo_pca)
  rbg_oppo_vars

  ## Check variables contributions to PCs:
  fviz_contrib(rgb_oppo_pca, choice = "var", axes = 1, top = 30)
  fviz_contrib(rgb_oppo_pca, choice = "var", axes = 2, top = 30)

## Create dataframe with scores:
scores_oppo = as.data.frame(rgb_oppo_pca$x)
str(scores_oppo)

indv_oppo = hsv_oppo[,c(1:4)]
str (indv_oppo)

scores_indv_oppo = cbind (indv_oppo,scores_oppo)
str(scores_indv_oppo)

## Plot results:
ggplot(scores_indv_oppo, aes(x=PC1, y=PC2, fill=MORPH, shape=MORPH)) +
  theme_bw() +
  geom_point(size = 6) +
  scale_shape_manual(name = "Season", values = c(21, 24)) +
  scale_fill_manual(name = "Morph", values=c("#AD7135","#cdb79e")) +
  xlab(paste0("PC1 (", round(perc_oppo_variance[1]*100,2),"%)")) +
  ylab(paste0("PC2 (", round(perc_oppo_variance[2]*100,2),"%)"))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size = 14)) +
  guides (shape = guide_legend (order=1), fill = guide_legend(order=2)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))


## Biplot of multivariate data:

  ## Add group columns:
  oppo_groups <- as.factor (paste (scores_indv_oppo$MORPH,scores_indv_oppo$SEASON, sep = " "))

  ## Plot Biplot:
  ggbiplot(rgb_oppo_pca, choices = 1:2, groups = oppo_groups, varname.size = 4, labels.size = 4, ellipse = TRUE,
          ellipse.alpha = 0) +
    theme_bw() +
    scale_color_manual(name = "Phenotypic group", values=c("#AD7135","#cdb79e")) +
    scale_fill_manual(name = "Phenotypic group", values=c("#AD7135","#cdb79e")) +
    scale_shape_manual(name = "Phenotypic group", values = c(21, 24)) +
    geom_point(aes(color= oppo_groups, fill=oppo_groups, shape=oppo_groups), size = 4) +
    xlab(paste0("PC1 (", round(perc_oppo_variance[1]*100,2),"%)")) +
    ylab(paste0("PC2 (", round(perc_oppo_variance[2]*100,2),"%)")) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size = 12)) +
    xlim(-2.2, 2.5)

## RUN KW TEST:
  
  ## Group the data by variable
  grouped.oppo <- hsv_oppo %>%
    gather(key = "variable", value = "value", H_h, UD_h, MD_h, LD_h, H_s, UD_s, MD_s, LD_s, H_v, UD_v, MD_v, LD_v) %>%
    group_by(variable)

  ## Do Kruskal-Wallis test
  grouped.oppo %>% kruskal_test(value ~ MORPH)



#### 6. BOXPLOTS OF EACH VARIABLE UNDER STUDY #####
#rm(list=ls())

#Read input: [if starting code here]
#m <- read.table("HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt", 
#                header=T, sep = "\t", na.strings = "n/a")
str (m)
table (m$MORPH, m$SEASON)


## Check columns for amount of missing data:
## No missing data left because rgb2hsv was applied to filtered dataset already
hsv_ind_total <- m[complete.cases(m), ]

str(hsv_ind_total)
summary (hsv_ind_total)


## Boxplot of colour datapoint by morph/season:


  ## Lower Dorsal:

  ggplot (hsv_ind_total, aes(x=MORPH, y=LD_h)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Lower Dorsum - Hue", ) +
    ylim (0.05,0.12) +
    theme(plot.title = element_text(hjust = 0.5))

  ggplot (hsv_ind_total, aes(x=MORPH, y=LD_s)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Lower Dorsum - Saturation", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))

  ggplot (hsv_ind_total, aes(x=MORPH, y=LD_v)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Lower Dorsum - Value", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))

  ## Middle Dorsal:

  ggplot (hsv_ind_total, aes(x=MORPH, y=MD_h)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Middle Dorsum - Hue", ) +
    ylim (0.05,0.12) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=MD_s)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Middle Dorsum - Saturation", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=MD_v)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Middle Dorsum - Value", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))
  

  ## Upper Dorsal:

  ggplot (hsv_ind_total, aes(x=MORPH, y=UD_h)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Upper Dorsum - Hue", ) +
    ylim (0.05,0.12) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=UD_s)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Upper Dorsum - Saturation", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=UD_v)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Upper Dorsum - Value", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))

  ## Head Dorsal:
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=H_h)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Head - Hue", ) +
    ylim (0.05,0.12) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=H_s)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Head - Saturation", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggplot (hsv_ind_total, aes(x=MORPH, y=H_v)) +
    theme_bw() +
    scale_fill_manual(name = "Season", values=c("#bebebe","#8b5a2b")) +
    geom_boxplot(aes(fill=fct_rev(SEASON)))+
    ggtitle ("Head - Value", ) +
    ylim (0,1) +
    theme(plot.title = element_text(hjust = 0.5))

