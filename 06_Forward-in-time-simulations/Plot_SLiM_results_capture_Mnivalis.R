## Plot SLiM results:
## Capture dataset analysis
## Last edited: 11-11-2025, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


#### Import data ####
data10 <- read.table("./ANC13_wMUT_neutral_presplit10perc_rescale5_final_400reps_STATS_phenoFreqs.txt", sep = "\t",
                   header = T)
str(data10)

data30 <- read.table("ANC13_wMUT_neutral_presplit30perc_rescale5_final_400reps_STATS_phenoFreqs.txt", sep = "\t",
                     header = T)
str(data30)

data50 <- read.table("ANC13_wMUT_neutral_presplit50perc_rescale5_final_400reps_STATS_phenoFreqs.txt", sep = "\t",
                     header = T)
str(data50)



#### Violin plots ####

## Import libraries:
library(dplyr) 
library(tidyverse)

## Get columns of interest
pheno10 <- data10[, c(13:16)]
str(pheno10)

pheno30 <- data30[, c(13:16)]
str(pheno30)

pheno50 <- data50[, c(13:16)]
str(pheno50)


## Convert data from wide format to long format and concatenate columns:
pheno10_a <- pivot_longer(pheno10, FreqBphenoP1, names_to = "FreqBeg", values_to = "FinalFreq")
pheno10_b <- pivot_longer(pheno10, FreqBphenoP2, names_to = "FreqBeg", values_to = "FinalFreq")

pheno30_a <- pivot_longer(pheno30, FreqBphenoP1, names_to = "FreqBeg", values_to = "FinalFreq")
pheno30_b <- pivot_longer(pheno30, FreqBphenoP2, names_to = "FreqBeg", values_to = "FinalFreq")

pheno50_a <- pivot_longer(pheno50, FreqBphenoP1, names_to = "FreqBeg", values_to = "FinalFreq")
pheno50_b <- pivot_longer(pheno50, FreqBphenoP2, names_to = "FreqBeg", values_to = "FinalFreq")


## Prepare data for concatenation:

  ## pre-split freq=10%:
  pheno10_a$Pheno <- "P1_Brown"
  pheno10_b$Pheno <- "P2_Brown"
  pheno10_a$Pop <- "P1"
  pheno10_b$Pop <- "P2"
  pheno10_a$FreqBeg <- "10perc"
  pheno10_b$FreqBeg <- "10perc"

    ## keep columns of interest:
    pheno10_a <- pheno10_a [, c(4:7)]
    pheno10_b <- pheno10_b [, c(4:7)]

  ## pre-split freq=30%:
  pheno30_a$Pheno <- "P1_Brown"
  pheno30_b$Pheno <- "P2_Brown"
  pheno30_a$Pop <- "P1"
  pheno30_b$Pop <- "P2"
  pheno30_a$FreqBeg <- "30perc"
  pheno30_b$FreqBeg <- "30perc"

    ## keep columns of interest:
    pheno30_a <- pheno30_a [, c(4:7)]
    pheno30_b <- pheno30_b [, c(4:7)]

  ## pre-split freq=50%:
  pheno50_a$Pheno <- "P1_Brown"
  pheno50_b$Pheno <- "P2_Brown"
  pheno50_a$Pop <- "P1"
  pheno50_b$Pop <- "P2"
  pheno50_a$FreqBeg <- "50perc"
  pheno50_b$FreqBeg <- "50perc"

    ## keep columns of interest:
    pheno50_a <- pheno50_a [, c(4:7)]
    pheno50_b <- pheno50_b [, c(4:7)]


## Concatenate data:
data_B_pheno <- as.data.frame(rbind(pheno10_a, pheno10_b,
                                    pheno30_a, pheno30_b,
                                    pheno50_a, pheno50_b)) 
str(data_B_pheno)


## Plot phenotype frequency in violin plots:

  ## Import library:
  library(ggbeeswarm)

  ## Violin plots scaled by width of density plot:
  ## Added horizontal lines with thresholds for simulations:
  ggplot(data_B_pheno, aes(x=FreqBeg, y=FinalFreq, fill=Pheno, color=Pheno)) +
    theme_bw () +
    geom_violin(scale = "width", width = 0.5, alpha = 0.1, position = position_dodge(width = 0.7)) +
    scale_color_manual(values=c("#DDAA33", "#004488")) +
    scale_fill_manual(values=c("#DDAA33", "#004488")) +
    geom_quasirandom(alpha = 0.3, width = 0.1, dodge.width = 0.7, varwidth = TRUE) +

    # --- Reset colour scale for the lines ---
    ggnewscale::new_scale_color() +
    
    # --- Plot polymorphism range limits ---  
    geom_hline(aes(yintercept = 0.05, color = "0.05 < f < 0.95"),   ## to plot hlines in legends, aes() is needed
              linetype = "dashed", , size = 1) +
    geom_hline(aes(yintercept = 0.20, color = "0.20 < f(vulgaris) < 0.80"),
              linetype = "dashed", size = 1) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey", size = 1) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey30", size = 1) +
    
    scale_color_manual(values=c("grey30", "grey")) +    ## manually set the color of hlines with aes()
    
    xlab("Pre-split f(vulgaris)") +
    ylab("Current f(vulgaris)") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.5))
  


#### Frequency of polymorphism maintenance ####

## Create new dataframe:
freqs10 <- data.frame(P1 = data10$FreqBphenoP1, P2 = data10$FreqBphenoP2)
freqs30 <- data.frame(P1 = data30$FreqBphenoP1, P2 = data30$FreqBphenoP2)
freqs50 <- data.frame(P1 = data50$FreqBphenoP1, P2 = data50$FreqBphenoP2)

## Add pre-split frequency:
freqs10$freqIni <- "10"
freqs30$freqIni <- "30"
freqs50$freqIni <- "50"

## Assess if final frequency of B/vulg within 0.2 to 0.8 range:
freqs10$inter <- freqs10$P1 >= 0.2 & freqs10$P1 <= 0.8 & freqs10$P2 >= 0.2 & freqs10$P2 <= 0.8
freqs30$inter <- freqs30$P1 >= 0.2 & freqs30$P1 <= 0.8 & freqs30$P2 >= 0.2 & freqs30$P2 <= 0.8
freqs50$inter <- freqs50$P1 >= 0.2 & freqs50$P1 <= 0.8 & freqs50$P2 >= 0.2 & freqs50$P2 <= 0.8

## Assess if final frequency of B/vulg within 0.05 to 0.95 range:
freqs10$large <- freqs10$P1 >= 0.05 & freqs10$P1 <= 0.95 & freqs10$P2 >= 0.05 & freqs10$P2 <= 0.95
freqs30$large <- freqs30$P1 >= 0.05 & freqs30$P1 <= 0.95 & freqs30$P2 >= 0.05 & freqs30$P2 <= 0.95
freqs50$large <- freqs50$P1 >= 0.05 & freqs50$P1 <= 0.95 & freqs50$P2 >= 0.05 & freqs50$P2 <= 0.95

## Concatenate dataset:
totalfreqs <- rbind (freqs10, freqs30, freqs50)

## Plot frequency of maintenance:

  ## Import library:
  library (patchwork)

  ## Plot results:
  p2 <- ggplot(totalfreqs, aes(x = factor(freqIni), fill = large)) +
    geom_bar(position = "fill") +
    ylab("Proportion of replicates") + xlab("Pre-split f(vulgaris)") +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_fill_manual(values = c("FALSE" = "#bebebe", "TRUE" = "#8b5a2b"),
                      name = "Maintenance of polymorphism at 0.05 - 0.95")

  p1 <- ggplot(totalfreqs, aes(x = factor(freqIni), fill = inter)) +
    geom_bar(position = "fill") +
    ylab("Proportion of replicates") + xlab("Pre-split f(vulgaris)") +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_fill_manual(values = c("FALSE" = "#bebebe", "TRUE" = "#8b5a2b"),
                      name = "Maintenance of polymorphism at 0.20 - 0.80")

  p1 / p2


