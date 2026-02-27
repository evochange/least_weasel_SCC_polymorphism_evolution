## MC1R haplotype divergence (Dxy):
## Capture dataset analysis
## Last edited: 10-05-2024, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())

## Prepare input files:
## Read in the file with sliding window estimates of Dxy:
  
  ## Genome-wide stats:
  stats_gwallHomoHG1HG2 <-read.csv("./capture_nivalis_allHomo_HG1HG2_2kfrags_regions.genomic.Dxy.popgenWind.csv",header=T)
  str(stats_gwallHomoHG1HG2)

  ## MC1R stats:
  stats_mc1rallHomoHG1HG2 <-read.csv("./capture_nivalis_allHomo_HG1HG2_mc1r_region.genomic.Dxy.popgenWind.csv",header=T)
  str(stats_mc1rallHomoHG1HG2)


## Import libraries:
library(dplyr) 
library(tidyverse)


## Draw plot - HG1_vulgaris vs HG2_vulgaris:

  ## calculate genome-wide mean:
  mean(stats_gwallHomoHG1HG2$dxy_homoB_HG1_homoB_HG2)  #0.0040362

  ## calculate genome-wide quantiles:
  quantile (stats_gwallHomoHG1HG2$dxy_homoB_HG1_homoB_HG2, probs = c(0.05, 0.50, 0.95, 0.99))
  #5%      50%      95%      99% 
  #0.001800 0.003900 0.007000 0.008401 

  ## plot data:
  ggplot(stats_mc1rallHomoHG1HG2, aes(x = mid, y = dxy_homoB_HG1_homoB_HG2))+
    geom_rect(xmin = 145924, xmax = 185924, ymin = 0, ymax = 0.020, fill = "grey90") +
    geom_point (aes(color = ifelse(dxy_homoB_HG1_homoB_HG2 > 0.007000 | dxy_homoB_HG1_homoB_HG2 < 0.001800,"black","darkorange"))) +
    scale_colour_manual(labels = c(">95%/<5% gw", "5%-95% gw"), values=c('darkorange', 'black')) +
    theme_light() +
    xlim (100000,250000) +
    ylim (0.00, 0.02) +
    xlab ("Scaffold position") + 
    ylab ("Dxy - HG1-vulg vs HG2-vulg") +
    labs (color = "Threshold") +
    geom_hline (yintercept = 0.007000, linetype = "dashed",color ="grey30", size = 1) +
    geom_hline (yintercept = 0.001800, linetype = "dashed",color ="grey30", size = 1) + 
    theme(
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )


## Draw plot - HG1_nivalis vs HG2_nivalis:

  ## calculate genome-wide mean:
  mean(stats_gwallHomoHG1HG2$dxy_homoW_HG1_homoW_HG2)  #0.0039433

  ## calculate genome-wide quantiles:
  quantile (stats_gwallHomoHG1HG2$dxy_homoW_HG1_homoW_HG2, probs = c(0.05, 0.50, 0.95, 0.99))
  #5%      50%      95%      99% 
  #0.001800 0.003700 0.007005 0.008400 

  ## plot data:
  ggplot(stats_mc1rallHomoHG1HG2, aes(x = mid, y = dxy_homoW_HG1_homoW_HG2))+
    geom_rect(xmin = 145924, xmax = 185924, ymin = 0, ymax = 0.020, fill = "grey90") +
    geom_point (aes(color = ifelse(dxy_homoW_HG1_homoW_HG2 > 0.007005 | dxy_homoW_HG1_homoW_HG2 < 0.001800,"black","darkorange"))) +
    scale_colour_manual(labels = c(">95%/<5% gw", "5%-95% gw"), values=c('darkorange', 'black')) +
    theme_light() +
    xlim (100000,250000) +
    ylim (0.00, 0.02) +
    xlab ("Scaffold position") + 
    ylab ("Dxy - HG1-niv vs HG2-niv") +
    labs (color = "Threshold") +
    geom_hline (yintercept = 0.007005, linetype = "dashed",color ="grey30", size = 1) +
    geom_hline (yintercept = 0.001800, linetype = "dashed",color ="grey30", size = 1) + 
    theme(
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )


## Draw plot - HG1_vulgaris vs HG1_nivalis:

  ## calculate genome-wide mean:
  mean(stats_gwallHomoHG1HG2$dxy_homoB_HG1_homoW_HG1)  #0.0037059

  ## calculate genome-wide quantiles:
  quantile (stats_gwallHomoHG1HG2$dxy_homoB_HG1_homoW_HG1, probs =  c(0.05, 0.50, 0.95, 0.99))
  #5%      50%      95%      99%  
  #0.001600 0.003500 0.006700 0.007702 

  ## plot data:
  ggplot(stats_mc1rallHomoHG1HG2, aes(x = mid, y = dxy_homoB_HG1_homoW_HG1))+
    geom_rect(xmin = 145924, xmax = 185924, ymin = 0, ymax = 0.020, fill = "grey90") +
    geom_point (aes(color = ifelse(dxy_homoB_HG1_homoW_HG1 > 0.006700 | dxy_homoB_HG1_homoW_HG1 < 0.001600,"black","darkorange"))) +
    scale_colour_manual(labels = c(">95%/<5% gw", "5%-95% gw"), values=c('darkorange', 'black')) +
    theme_light() +
    xlim (100000,250000) +
    ylim (0.00, 0.02) +
    xlab ("Scaffold position") + 
    ylab ("Dxy - HG1-vulg vs HG1-niv") +
    labs (color = "Threshold") +
    geom_hline (yintercept = 0.006700, linetype = "dashed",color ="grey30", size = 1) +
    geom_hline (yintercept = 0.001600, linetype = "dashed",color ="grey30", size = 1) + 
    theme(
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )


## Draw plot - HG2_vulgaris vs HG2_nivalis:

  ## calculate genome-wide mean:
  mean(stats_gwallHomoHG1HG2$dxy_homoB_HG2_homoW_HG2)  #0.0033346

  ## calculate genome-wide quantiles:
  quantile (stats_gwallHomoHG1HG2$dxy_homoB_HG2_homoW_HG2, probs =  c(0.05, 0.50, 0.95, 0.99))
  #5%      50%      95%      99%  
  #0.001295 0.003100 0.006205 0.007600 

  ## plot data:
  ggplot(stats_mc1rallHomoHG1HG2, aes(x = mid, y = dxy_homoB_HG2_homoW_HG2))+
    geom_rect(xmin = 145924, xmax = 185924, ymin = 0, ymax = 0.020, fill = "grey90") +
    geom_point (aes(color = ifelse(dxy_homoB_HG2_homoW_HG2 > 0.006205 | dxy_homoB_HG2_homoW_HG2 < 0.001295,"black","darkorange"))) +
    scale_colour_manual(labels = c(">95%/<5% gw", "5%-95% gw"), values=c('darkorange', 'black')) +
    theme_light() +
    xlim (100000,250000) +
    ylim (0.00, 0.02) +
    xlab ("Scaffold position") + 
    ylab ("Dxy - HG2-vulg vs HG2-niv") +
    labs (color = "Threshold") +
    geom_hline (yintercept = 0.006205, linetype = "dashed",color ="grey30", size = 1) +
    geom_hline (yintercept = 0.001295, linetype = "dashed",color ="grey30", size = 1) + 
    theme(
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )

  