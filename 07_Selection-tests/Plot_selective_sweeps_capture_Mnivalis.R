## Plot selective sweep analysis:
## Capture dataset analysis
## Last edited: 15-03-2024, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


## Import libraries:
require (ggplot2)
require (scales)


#### Dataset: vulgaris morph ####

## Read data:
sweed_vulg <- read.table("./SweeD_Report.capture_snpsVULGARIS_scaf337.strictPol", header = T)

## Plot results:
ggplot(data = sweed_vulg, mapping = aes(x=Position, y=Likelihood)) +
  theme_classic() +
  geom_col (color="black", fill="black",width = 10000) +
  scale_y_continuous(ylab("CLR"), expand = c(0.005, 0.005),breaks =c(0, 250, 500), limits = c(0,500)) +
  scale_x_continuous(xlab ("Position"), breaks =c(0,250000,500000), limits = c(0,510000), expand = c(0.005, 0.005)) +
  geom_segment(aes(x=0, xend=500000, y=70.81, yend=70.81), colour="red",size=1, lty="dashed") +
  ggtitle ("vulgaris morph") +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=14)
  )

  
#### Dataset: nivalis morph ####

## Read data:
sweed_niv <- read.table("./SweeD_Report.capture_snpsNIVALIS_scaf337.strictPol", header = T)

## Plot results:
ggplot(data = sweed_niv, mapping = aes(x=Position, y=Likelihood)) +
  theme_classic() +
  geom_col (color="black", fill="black",width = 10000) +
  scale_y_continuous(ylab("CLR"), expand = c(0.005, 0.005),breaks =c(0, 250, 500), limits = c(0,500)) +
  scale_x_continuous(xlab ("Position"), breaks =c(0,250000,500000), limits = c(0,510000), expand = c(0.005, 0.005)) +
  geom_segment(aes(x=0, xend=500000, y=70.81, yend=70.81), colour="red",size=1, lty="dashed") +
  ggtitle ("nivalis morph") +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=14)
  )


#### Dataset: homozygous HG1-vulg ####

## Read data:
sweed_HG1_B <- read.table("./SweeD_Report.capture_homoBrown_HG1_scaf337.strictPol", header = T)

## Plot results:
ggplot(data = sweed_HG1_B, mapping = aes(x=Position, y=Likelihood)) +
  theme_classic() +
  geom_col (color="black", fill="black",width = 10000) +
  scale_y_continuous(ylab("CLR"), expand = c(0.005, 0.005),breaks =c(0, 250, 500), limits = c(0,500)) +
  scale_x_continuous(xlab ("Position"), breaks =c(0,250000,500000), limits = c(0,510000), expand = c(0.005, 0.005)) +
  geom_segment(aes(x=0, xend=500000, y=8.17, yend=8.17), colour="red",size=1, lty="dashed") +
  ggtitle ("homozygous HG1-vulg") +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=14)
  )


#### Dataset: homozygous HG1-niv ####

## Read data:
sweed_HG1_W <- read.table("./SweeD_Report.capture_homoWhite_HG1_scaf337.strictPol", header = T)

## Plot results:
ggplot(data = sweed_HG1_W, mapping = aes(x=Position, y=Likelihood)) +
  theme_classic() +
  geom_col (color="black", fill="black",width = 10000) +
  scale_y_continuous(ylab("CLR"), expand = c(0.005, 0.005),breaks =c(0, 250, 500), limits = c(0,500)) +
  scale_x_continuous(xlab ("Position"), breaks =c(0,250000,500000), limits = c(0,510000), expand = c(0.005, 0.005)) +
  geom_segment(aes(x=0, xend=500000, y=8.17, yend=8.17), colour="red",size=1, lty="dashed") +
  ggtitle ("homozygous HG1-niv") +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=14)
  )


#### Dataset: homozygous HG2-vulg ####

## Read data:
sweed_HG2_B <- read.table("./SweeD_Report.capture_homoBrown_HG2_scaf337.strictPol", header = T)

## Plot results:
ggplot(data = sweed_HG2_B, mapping = aes(x=Position, y=Likelihood)) +
  theme_classic() +
  geom_col (color="black", fill="black",width = 10000) +
  scale_y_continuous(ylab("CLR"), expand = c(0.005, 0.005),breaks =c(0, 250, 500), limits = c(0,500)) +
  scale_x_continuous(xlab ("Position"), breaks =c(0,250000,500000), limits = c(0,510000), expand = c(0.005, 0.005)) +
  geom_segment(aes(x=0, xend=500000, y=8.17, yend=8.17), colour="red",size=1, lty="dashed") +
  ggtitle ("homozygous HG2-vulg") +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=14)
  )


#### Dataset: homozygous HG2-niv ####

## Read data:
sweed_HG2_W <- read.table("./SweeD_Report.capture_homoWhite_HG2_scaf337.strictPol", header = T)

## Plot results:
ggplot(data = sweed_HG2_W, mapping = aes(x=Position, y=Likelihood)) +
  theme_classic() +
  geom_col (color="black", fill="black",width = 10000) +
  scale_y_continuous(ylab("CLR"), expand = c(0.005, 0.005),breaks =c(0, 250, 500), limits = c(0,500)) +
  scale_x_continuous(xlab ("Position"), breaks =c(0,250000,500000), limits = c(0,510000), expand = c(0.005, 0.005)) +
  geom_segment(aes(x=0, xend=500000, y=8.17, yend=8.17), colour="red",size=1, lty="dashed") +
  ggtitle ("homozygous HG2-niv") +
  theme(
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    plot.title = element_text(size=14)
  )

