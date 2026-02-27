## Plot MapDamage results:
## Capture dataset analysis
## Last edited: 26-09-2023, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


## Import data to plot:
data <- read.table("./capture_Mnivalis_other_mustelids_mapDamage_M_SD.txt", header = T,
                   sep = "\t")
head (data)


## Plot results:

  ## Import library:
  library (ggplot2)

  ## Draw plot:
  ggplot(data, aes(x=Position, y=Mean, ymin=(Mean-SD), ymax=(Mean+SD), 
                  fill=Dataset, linetype=Dataset)) +
    theme_bw() + 
    geom_line (size = 1) +
    geom_ribbon (alpha = 0.5) +
    scale_fill_manual(values = c("#66CCEE", "#AA3377", "#CCBB44")) +
    ylab("C-to-T Substitution Frequency") +
    xlab ("Distance from read end (bp)") + 
    ylim (0, 0.03)
  