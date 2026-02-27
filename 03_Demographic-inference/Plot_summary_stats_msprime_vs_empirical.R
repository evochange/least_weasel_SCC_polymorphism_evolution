## PLOT msprime simulation results
## Demography model-based stats
## VS EMPIRICAL DATA
## Last edited: 30-12-2024, Inês Miranda


#set working directory:
setwd("C:/path/to/folder/")
rm(list = ls())


#### 1. Empirical Data ####

##### 1.1. Tajima's D data - Estimates per population #####

## Read in the files with sliding window estimates of TajD:
## Genome-wide stats:
stats_gwNorth<-read.csv("./capture_europe_North_2k_regions.genomic.TajD.popFreq.csv", header=T)
str(stats_gwNorth)
stats_gwCentral<-read.csv("./capture_europe_Central_2k_regions.genomic.TajD.popFreq.csv", header=T)
str(stats_gwCentral)

## Extract TajD values from table:

  ## import libraries:
  library(dplyr) 
  library(tidyverse)

  ## Extract stats:
  tajWind_gwNorth <- stats_gwNorth %>% select(contains("TajD_"))
  str(tajWind_gwNorth)
  tajWind_gwCentral <- stats_gwCentral %>% select(contains("TajD_"))
  str(tajWind_gwCentral)

  ## Rename column names:
  colnames(tajWind_gwNorth)[1] <-  "TajD_North_empirical"
  colnames(tajWind_gwCentral)[1] <-  "TajD_Central_empirical"


##### 1.2 Fst data - Pairwise estimations #####

## Read in the file with sliding window estimates of Fst:
stats_gw13<-read.csv("./capture_europe_NC_2k_regions_genomic.Fst.popgenWind.csv",header=T)
str(stats_gw13)

## Extract Fst values from table:

  ## Extract stats:
  fst_gw13 <- stats_gw13 %>% select(contains("Fst_"))
  str(fst_gw13)

  ## Rename column name:
  colnames(fst_gw13)[1] <-  "Fst_NC_empirical"
  str (fst_gw13)

  
#### 2. MSPrime Simulations Data ####

## Read in file with msprime simulated values for TajD and Fst:
sims_ANC13 <- read.table("./20240223_msprime_sims_NorthCentral_gadmaBest_Fst_TajD_1k_reps.txt",
                             header = T, sep = "\t")
str(sims_ANC13)
    
    
#### 3. Plot simulations vs empirical data - North/Central model ####

## 3.1. TajD data - North:

## Combined empirical and simulated data in a single dataframe:
taj_North_13 <- as.data.frame(cbind(tajWind_gwNorth$TajD_North_empirical, 
                                      sims_ANC13$TajD_North)) 
colnames (taj_North_13) <- c("North_empirical", "North_sims13")
str(taj_North_13)

taj_North_emp <- pivot_longer(taj_North_13[1], North_empirical, names_to = "group", values_to = "pi")
taj_North_sims13 <- pivot_longer(taj_North_13[2], North_sims13, names_to = "group", values_to = "pi")

taj_North_13_total <- as.data.frame(rbind (taj_North_emp, taj_North_sims13))
str(taj_North_13_total)

    ## Draw density plots - North:
    ggplot(data = taj_North_13_total, 
             aes(x = pi, group = group, color = group, linetype = group)) +
      theme_bw() +
      geom_density(size = 1) +
      ylab ("Density") +
      xlab ("North - Taj D")
    
    
## 3.2. TajD data - Central:

## Combined empirical and simulated data in a single dataframe:
taj_Central_13 <- as.data.frame(cbind(tajWind_gwCentral$TajD_Central_empirical, 
                                      sims_ANC13$TajD_Central)) 
colnames (taj_Central_13) <- c("Central_empirical", "Central_sims13")
str(taj_Central_13)
    
taj_Central_emp <- pivot_longer(taj_Central_13[1], Central_empirical, names_to = "group", values_to = "pi")
taj_Central_sims13 <- pivot_longer(taj_Central_13[2], Central_sims13, names_to = "group", values_to = "pi")

taj_Central_sims13_total <- as.data.frame(rbind (taj_Central_emp,taj_Central_sims13))
str(taj_Central_sims13)
    
    ## Draw density plots:
    ggplot(data = taj_Central_sims13_total,
           aes(x = pi, group = group, color = group, linetype = group)) +
      theme_bw() +
      geom_density(size = 1) +
      ylab ("Density") +
      xlab ("Central - Taj D")


## 3.3. Fst data - North - Central:

## Combined empirical and simulated data in a single dataframe:
fst_NC_13 <- as.data.frame(cbind(fst_gw13$Fst_NC_empirical, 
                                 sims_ANC13$Fst_NorthCentral)) 
colnames (fst_NC_13) <- c("NC_empirical", "NC_sims13")
str(fst_NC_13)

fst_NC_emp <- pivot_longer(fst_NC_13[1], NC_empirical, names_to = "group", values_to = "pi")
fst_NC_sims13 <- pivot_longer(fst_NC_13[2], NC_sims13, names_to = "group", values_to = "pi")

fst_NC_sims13 <- as.data.frame(rbind (fst_NC_emp,fst_NC_sims13))
str(fst_NC_sims13)

    ## Draw density plot:
    ggplot(data = fst_NC_sims13, 
           aes(x = pi, group = group, color = group, linetype = group)) +
      theme_bw() +
      geom_density(size = 1) +
      ylab ("Density") +
      xlab ("FST - North-Central")
