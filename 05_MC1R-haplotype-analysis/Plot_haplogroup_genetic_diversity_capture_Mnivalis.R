## MC1R genetic diversity (per haplogroup):
## Capture dataset analysis
## Last edited: 22-05-2024, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


#### Pi Boxplots ####


## Prepare input files:
## Read in files with sliding window estimates of Pi:

  ## Genome-wide stats:
  stats_gwbrownHG1<-read.csv("./capture_nivalis_homoBrown_HG1_2kfrags_regions.genomic.popgenWind.csv",header=T)
  str(stats_gwbrownHG1)
  stats_gwwhiteHG1<-read.csv("./capture_nivalis_homoWhite_HG1_2kfrags_regions.genomic.popgenWind.csv",header=T)
  str(stats_gwwhiteHG1)
  stats_gwbrownHG2<-read.csv("./capture_nivalis_homoBrown_HG2_2kfrags_regions.genomic.popgenWind.csv",header=T)
  str(stats_gwbrownHG2)
  stats_gwwhiteHG2<-read.csv("./capture_nivalis_homoWhite_HG2_2kfrags_regions.genomic.popgenWind.csv",header=T)
  str(stats_gwwhiteHG2)
  stats_gwEuroALL<-read.csv("./capture_europeALL_2kfrags_regions.genomic.popgenWind.csv",header=T)
  str(stats_gwEuroALL)

  ## MC1R stats:
  stats_mc1rbrownHG1<-read.csv("./capture_nivalis_homoBrown_HG1_mc1r_40kb_region.genomic.popgenWind.csv",header=T)
  str(stats_mc1rbrownHG1)
  stats_mc1rwhiteHG1<-read.csv("./capture_nivalis_homoWhite_HG1_mc1r_40kb_region.genomic.popgenWind.csv",header=T)
  str(stats_mc1rwhiteHG1)
  stats_mc1rbrownHG2<-read.csv("./capture_nivalis_homoBrown_HG2_mc1r_40kb_region.genomic.popgenWind.csv",header=T)
  str(stats_mc1rbrownHG2)
  stats_mc1rwhiteHG2<-read.csv("./capture_nivalis_homoWhite_HG2_mc1r_40kb_region.genomic.popgenWind.csv",header=T)
  str(stats_mc1rwhiteHG2)
  stats_mc1rEuroALL<-read.csv("./capture_europeALL_mc1r_40kb_regkion.genomic.popgenWind.csv",header=T)
  str(stats_mc1rEuroALL)


## Extract pi values from table:

  ## Import libraries:
  library(dplyr) 
  library(tidyverse)

  ## GW stats:
  piWind_gwbrownHG1 <- stats_gwbrownHG1 %>% select(contains("pi_"))
  str(piWind_gwbrownHG1)
  piWind_gwwhiteHG1 <- stats_gwwhiteHG1 %>% select(contains("pi_"))
  str(piWind_gwwhiteHG1)
  piWind_gwbrownHG2 <- stats_gwbrownHG2 %>% select(contains("pi_"))
  str(piWind_gwbrownHG2)
  piWind_gwwhiteHG2 <- stats_gwwhiteHG2 %>% select(contains("pi_"))
  str(piWind_gwwhiteHG2)
  piWind_gwEuroALL <- stats_gwEuroALL %>% select(contains("pi_"))
  str(piWind_gwEuroALL)

  ## MC1R stats:
  piWind_mc1rbrownHG1 <- stats_mc1rbrownHG1 %>% select(contains("pi_"))
  str(piWind_mc1rbrownHG1)
  piWind_mc1rwhiteHG1 <- stats_mc1rwhiteHG1 %>% select(contains("pi_"))
  str(piWind_mc1rwhiteHG1)
  piWind_mc1rbrownHG2 <- stats_mc1rbrownHG2 %>% select(contains("pi_"))
  str(piWind_mc1rbrownHG2)
  piWind_mc1rwhiteHG2 <- stats_mc1rwhiteHG2 %>% select(contains("pi_"))
  str(piWind_mc1rwhiteHG2)
  piWind_mc1rEuroALL <- stats_mc1rEuroALL %>% select(contains("pi_"))
  str(piWind_mc1rEuroALL)

  ## Rename column names:
  colnames(piWind_gwbrownHG1)[1] <-  "homoB_HG1"
  colnames(piWind_gwwhiteHG1)[1] <-  "homoW_HG1"
  colnames(piWind_gwbrownHG2)[1] <-  "homoB_HG2"
  colnames(piWind_gwwhiteHG2)[1] <-  "homoW_HG2"
  colnames(piWind_gwEuroALL)[1] <-  "EuropeALL"
  
  colnames(piWind_mc1rbrownHG1)[1] <-  "homoB_HG1"
  colnames(piWind_mc1rwhiteHG1)[1] <-  "homoW_HG1"
  colnames(piWind_mc1rbrownHG2)[1] <-  "homoB_HG2"
  colnames(piWind_mc1rwhiteHG2)[1] <-  "homoW_HG2"
  colnames(piWind_mc1rEuroALL)[1] <-  "EuropeALL"
  

## Convert data from wide format to long format:

## GW stats:
piWind_gwbrownHG1_2 <- pivot_longer(piWind_gwbrownHG1, homoB_HG1, names_to = "group", values_to = "Pi")
piWind_gwwhiteHG1_2 <- pivot_longer(piWind_gwwhiteHG1, homoW_HG1, names_to = "group", values_to = "Pi")
piWind_gwbrownHG2_2 <- pivot_longer(piWind_gwbrownHG2, homoB_HG2, names_to = "group", values_to = "Pi")
piWind_gwwhiteHG2_2 <- pivot_longer(piWind_gwwhiteHG2, homoW_HG2, names_to = "group", values_to = "Pi")
piWind_gwEuroALL_2 <- pivot_longer(piWind_gwEuroALL, EuropeALL, names_to = "group", values_to = "Pi")

## MC1R stats:
piWind_mc1rbrownHG1_2 <- pivot_longer(piWind_mc1rbrownHG1, homoB_HG1, names_to = "group", values_to = "Pi")
piWind_mc1rwhiteHG1_2 <- pivot_longer(piWind_mc1rwhiteHG1, homoW_HG1, names_to = "group", values_to = "Pi")
piWind_mc1rbrownHG2_2 <- pivot_longer(piWind_mc1rbrownHG2, homoB_HG2, names_to = "group", values_to = "Pi")
piWind_mc1rwhiteHG2_2 <- pivot_longer(piWind_mc1rwhiteHG2, homoW_HG2, names_to = "group", values_to = "Pi")
piWind_mc1rEuroALL_2 <- pivot_longer(piWind_mc1rEuroALL, EuropeALL, names_to = "group", values_to = "Pi")


## Combine into single dataframe per dataset (GW vs MC1r):

  ## GW stats:
  piWind_gwTotal <- as.data.frame(rbind(piWind_gwbrownHG1_2, piWind_gwwhiteHG1_2,
                                        piWind_gwbrownHG2_2, piWind_gwwhiteHG2_2,
                                        piWind_gwEuroALL_2)) 

    ## Add info from data type (genome-wide vs mc1r) in new column:
    piWind_gwTotal$data <- "genome-wide"
    str (piWind_gwTotal)

  ## MC1R stats:
  piWind_mc1rTotal <- as.data.frame(rbind(piWind_mc1rbrownHG1_2, piWind_mc1rwhiteHG1_2,
                                          piWind_mc1rbrownHG2_2, piWind_mc1rwhiteHG2_2,
                                          piWind_mc1rEuroALL_2))

    ## Add info from data type (genome-wide vs mc1r) in new column:
    piWind_mc1rTotal$data <- "mc1r"
    str (piWind_mc1rTotal)


## Stack the two dataframes (GW and MC1r) together:
piWind_Full <- as.data.frame(rbind (piWind_gwTotal, piWind_mc1rTotal))
str(piWind_Full)

## Draw boxplots:
ggplot(piWind_Full, aes(x = Pi, y = group, fill = data))+
  geom_boxplot () +
  theme_light() +
  ylim ("homoB_HG1", "homoW_HG1", "homoB_HG2", "homoW_HG2", "EuropeALL")+
  xlim (0, 0.012) +
  coord_flip () +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )


##### Run significance tests #####

## Import library:
library(rstatix)

## homoB_HG1:
piWind_brownHG1 <- piWind_Full [piWind_Full$group == "homoB_HG1",]
str(piWind_brownHG1)

  ## Do Kruskal-Wallis test
  boxplot (Pi ~ data, data = piWind_brownHG1)
  kruskal.test (Pi ~ data, data = piWind_brownHG1)
  piWind_brownHG1 %>% kruskal_test(Pi ~ data)

## homoB_HG2:
piWind_brownHG2 <- piWind_Full [piWind_Full$group == "homoB_HG2",]
str(piWind_brownHG2)

  ## Do Kruskal-Wallis test
  boxplot (Pi ~ data, data = piWind_brownHG2)
  kruskal.test (Pi ~ data, data = piWind_brownHG2)
  piWind_brownHG2 %>% kruskal_test(Pi ~ data) 

## white_HG1:
piWind_whiteHG1 <- piWind_Full [piWind_Full$group == "homoW_HG1",]
str(piWind_whiteHG1)

  ## Do Kruskal-Wallis test
  boxplot (Pi ~ data, data = piWind_whiteHG1)
  kruskal.test (Pi ~ data, data = piWind_whiteHG1)
  piWind_whiteHG1 %>% kruskal_test(Pi ~ data) 

## white_HG2:
piWind_whiteHG2 <- piWind_Full [piWind_Full$group == "homoW_HG2",]
str(piWind_whiteHG2)

  ## Do Kruskal-Wallis test
  boxplot (Pi ~ data, data = piWind_whiteHG2)
  kruskal.test (Pi ~ data, data = piWind_whiteHG2)
  piWind_whiteHG2 %>% kruskal_test(Pi ~ data) 

## EuropeALL:
piWind_ALL <- piWind_Full [piWind_Full$group == "EuropeALL",]
str(piWind_ALL)

  ## Do Kruskal-Wallis test
  boxplot (Pi ~ data, data = piWind_ALL)
  kruskal.test (Pi ~ data, data = piWind_ALL)
  piWind_ALL %>% kruskal_test(Pi ~ data) 



#### Tajima's D Boxplots ####

## Clean environment:
rm(list = ls())

## Prepare input files:
## Read in files with sliding window estimates of Tajima's D:

  ## Genome-wide stats:
  stats_gwbrownHG1<-read.csv("./capture_nivalis_homoBrown_HG1_2kfrags_regions.genomic.popFreq.csv",header=T)
  str(stats_gwbrownHG1)
  stats_gwwhiteHG1<-read.csv("./capture_nivalis_homoWhite_HG1_2kfrags_regions.genomic.popFreq.csv",header=T)
  str(stats_gwwhiteHG1)
  stats_gwbrownHG2<-read.csv("./capture_nivalis_homoBrown_HG2_2kfrags_regions.genomic.popFreq.csv",header=T)
  str(stats_gwbrownHG2)
  stats_gwwhiteHG2<-read.csv("./capture_nivalis_homoWhite_HG2_2kfrags_regions.genomic.popFreq.csv",header=T)
  str(stats_gwwhiteHG2)
  stats_gwEuroALL<-read.csv("./capture_europeALL_2kfrags_regions.genomic.popFreq.csv",header=T)
  str(stats_gwEuroALL)

  ## MC1R stats:
  stats_mc1rbrownHG1<-read.csv("./capture_nivalis_homoBrown_HG1_mc1r_40kb_region.genomic.popFreq.csv",header=T)
  str(stats_mc1rbrownHG1)
  stats_mc1rwhiteHG1<-read.csv("./capture_nivalis_homoWhite_HG1_mc1r_40kb_region.genomic.popFreq.csv",header=T)
  str(stats_mc1rwhiteHG1)
  stats_mc1rbrownHG2<-read.csv("./capture_nivalis_homoBrown_HG2_mc1r_40kb_region.genomic.popFreq.csv",header=T)
  str(stats_mc1rbrownHG2)
  stats_mc1rwhiteHG2<-read.csv("./capture_nivalis_homoWhite_HG2_mc1r_40kb_region.genomic.popFreq.csv",header=T)
  str(stats_mc1rwhiteHG2)
  stats_mc1rEuroALL<-read.csv("./capture_europeALL_mc1r_40kb_region.genomic.popFreq.csv",header=T)
  str(stats_mc1rEuroALL)


## Extract TajD values from table:
  
  ## GW stats:
  tjaD_gwbrownHG1 <- stats_gwbrownHG1 %>% select(contains("TajD_"))
  str(tjaD_gwbrownHG1)
  tajD_gwwhiteHG1 <- stats_gwwhiteHG1 %>% select(contains("TajD_"))
  str(tajD_gwwhiteHG1)
  tajD_gwbrownHG2 <- stats_gwbrownHG2 %>% select(contains("TajD_"))
  str(tajD_gwbrownHG2)
  tajD_gwwhiteHG2 <- stats_gwwhiteHG2 %>% select(contains("TajD_"))
  str(tajD_gwwhiteHG2)
  tajD_gwEuroALL <- stats_gwEuroALL %>% select(contains("TajD_"))
  str(tajD_gwEuroALL)

  ## MC1R stats:
  tajD_mc1rbrownHG1 <- stats_mc1rbrownHG1 %>% select(contains("TajD_"))
  str(tajD_mc1rbrownHG1)
  tajD_mc1rwhiteHG1 <- stats_mc1rwhiteHG1 %>% select(contains("TajD_"))
  str(tajD_mc1rwhiteHG1)
  tajD_mc1rbrownHG2 <- stats_mc1rbrownHG2 %>% select(contains("TajD_"))
  str(tajD_mc1rbrownHG2)
  tajD_mc1rwhiteHG2 <- stats_mc1rwhiteHG2 %>% select(contains("TajD_"))
  str(tajD_mc1rwhiteHG2)
  tajD_mc1rEuroALL <- stats_mc1rEuroALL %>% select(contains("TajD_"))
  str(tajD_mc1rEuroALL)

  ## Rename column names:
  colnames(tjaD_gwbrownHG1)[1] <-  "homoB_HG1"
  colnames(tajD_gwwhiteHG1)[1] <-  "homoW_HG1"
  colnames(tajD_gwbrownHG2)[1] <-  "homoB_HG2"
  colnames(tajD_gwwhiteHG2)[1] <-  "homoW_HG2"
  colnames(tajD_gwEuroALL)[1] <-  "EuropeALL"
  
  colnames(tajD_mc1rbrownHG1)[1] <-  "homoB_HG1"
  colnames(tajD_mc1rwhiteHG1)[1] <-  "homoW_HG1"
  colnames(tajD_mc1rbrownHG2)[1] <-  "homoB_HG2"
  colnames(tajD_mc1rwhiteHG2)[1] <-  "homoW_HG2"
  colnames(tajD_mc1rEuroALL)[1] <-  "EuropeALL"


## Convert data from wide format to long format:

## GW stats:
tjaD_gwbrownHG1_2 <- pivot_longer(tjaD_gwbrownHG1, homoB_HG1, names_to = "group", values_to = "TajD")
tajD_gwwhiteHG1_2 <- pivot_longer(tajD_gwwhiteHG1, homoW_HG1, names_to = "group", values_to = "TajD")
tajD_gwbrownHG2_2 <- pivot_longer(tajD_gwbrownHG2, homoB_HG2, names_to = "group", values_to = "TajD")
tajD_gwwhiteHG2_2 <- pivot_longer(tajD_gwwhiteHG2, homoW_HG2, names_to = "group", values_to = "TajD")
tajD_gwEuroALL_2 <- pivot_longer(tajD_gwEuroALL, EuropeALL, names_to = "group", values_to = "TajD")

## MC1R stats:
tajD_mc1rbrownHG1_2 <- pivot_longer(tajD_mc1rbrownHG1, homoB_HG1, names_to = "group", values_to = "TajD")
tajD_mc1rwhiteHG1_2 <- pivot_longer(tajD_mc1rwhiteHG1, homoW_HG1, names_to = "group", values_to = "TajD")
tajD_mc1rbrownHG2_2 <- pivot_longer(tajD_mc1rbrownHG2, homoB_HG2, names_to = "group", values_to = "TajD")
tajD_mc1rwhiteHG2_2 <- pivot_longer(tajD_mc1rwhiteHG2, homoW_HG2, names_to = "group", values_to = "TajD")
tajD_mc1rEuroALL_2 <- pivot_longer(tajD_mc1rEuroALL, EuropeALL, names_to = "group", values_to = "TajD")


## Combine into single dataframe per dataset (GW vs MC1r):

  ## GW stats:
  tajD_gwTotal <- as.data.frame(rbind(tjaD_gwbrownHG1_2, tajD_gwwhiteHG1_2,
                                      tajD_gwbrownHG2_2, tajD_gwwhiteHG2_2,
                                      tajD_gwEuroALL_2)) 

    ## Add info from data type (genome-wide vs mc1r) in new column:
    tajD_gwTotal$data <- "genome-wide"
    str (tajD_gwTotal)

  ## MC1R stats:
  tajD_mc1rTotal <- as.data.frame(rbind(tajD_mc1rbrownHG1_2, tajD_mc1rwhiteHG1_2,
                                        tajD_mc1rbrownHG2_2, tajD_mc1rwhiteHG2_2,
                                        tajD_mc1rEuroALL_2))

    ## Add info from data type (genome-wide vs mc1r) in new column:
    tajD_mc1rTotal$data <- "mc1r"
    str (tajD_mc1rTotal)


## Stack the two dataframes (GW and MC1r) together:
tajD_Full <- as.data.frame(rbind (tajD_gwTotal, tajD_mc1rTotal))
str(tajD_Full)


## Draw boxplots:
ggplot(tajD_Full, aes(x = TajD, y = group, fill = data))+
  geom_boxplot () +
  theme_light() +
  xlim (-2.5, 2.5) +
  ylim ("homoB_HG1", "homoW_HG1", "homoB_HG2", "homoW_HG2", "EuropeALL")+
  coord_flip () +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )


##### Run significance tests #####

## homoB_HG1:
tajD_brownHG1 <- tajD_Full [tajD_Full$group == "homoB_HG1",]
str(tajD_brownHG1)

  ## Do Kruskal-Wallis test
  boxplot (TajD ~ data, data = tajD_brownHG1)
  kruskal.test (TajD ~ data, data = tajD_brownHG1)
  tajD_brownHG1 %>% kruskal_test(TajD ~ data)

## homoB_HG2:
tajD_brownHG2 <- tajD_Full [tajD_Full$group == "homoB_HG2",]
str(tajD_brownHG2)

  ## Do Kruskal-Wallis test
  boxplot (TajD ~ data, data = tajD_brownHG2)
  kruskal.test (TajD ~ data, data = tajD_brownHG2)
  tajD_brownHG2 %>% kruskal_test(TajD ~ data)

## homoW_HG1:
tajD_whiteHG1 <- tajD_Full [tajD_Full$group == "homoW_HG1",]
str(tajD_whiteHG1)

  ## Do Kruskal-Wallis test
  boxplot (TajD ~ data, data = tajD_whiteHG1)
  kruskal.test (TajD ~ data, data = tajD_whiteHG1)
  tajD_whiteHG1 %>% kruskal_test(TajD ~ data)

## homoW_HG2:
tajD_whiteHG2 <- tajD_Full [tajD_Full$group == "homoW_HG2",]
str(tajD_whiteHG2)

  ## Do Kruskal-Wallis test
  boxplot (TajD ~ data, data = tajD_whiteHG2)
  kruskal.test (TajD ~ data, data = tajD_whiteHG2)
  tajD_whiteHG2 %>% kruskal_test(TajD ~ data)

## EuropeALL:
tajD_ALL <- tajD_Full [tajD_Full$group == "EuropeALL",]
str(tajD_ALL)

  ## Do Kruskal-Wallis test
  boxplot (TajD ~ data, data = tajD_ALL)
  kruskal.test (TajD ~ data, data = tajD_ALL)
  tajD_ALL %>% kruskal_test(TajD ~ data) 

