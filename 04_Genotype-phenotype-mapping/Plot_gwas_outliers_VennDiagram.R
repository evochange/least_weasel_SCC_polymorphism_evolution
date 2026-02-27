## Plot Case Control - Venn diagram - overlaps:
## Capture dataset analysis
## Last edited: 24-10-2023, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


## Import outlier tables - Case control:

ccMOPRHS <- read.table ("./capture_gwasMORPHS_capture_Mnivalis.flt1.CC.BONF.valid.significantP05", 
                      sep = "\t", header = F)
colnames(ccMOPRHS) <- c("SNP_ID", "Pval")
str (ccMOPRHS)

fetMORPHS <- read.table ("./capture_gwasMORPHS_capture_Mnivalis.flt1.FET.significantP05", 
                       sep = "\t", header = F)
colnames(fetMORPHS) <- c("SNP_ID", "Pval")
str (fetMORPHS)

trendMORPHS <- read.table ("./capture_gwasMORPHS_capture_Mnivalis.flt1.CA.TREND.significantP05", 
                         sep = "\t", header = F)
colnames(trendMORPHS) <- c("SNP_ID", "Pval")
str (trendMORPHS)


## Create list with case-control outlier datasets to compare:
ccOutliers <-  list(CC = ccMOPRHS$SNP_ID,
                    FET = fetMORPHS$SNP_ID,
                    TREND = trendMORPHS$SNP_ID)

## Import packages and plot Venn:
library("ggVennDiagram")
library("ggplot2")


ggVennDiagram(ccOutliers, category.names = c("CC","FET","TREND"), 
              label_alpha = 0) +
  scale_color_brewer(palette = "Blues") +
  ggplot2::scale_fill_gradient(low="aliceblue",high = "cadetblue4") +
  ggtitle("Coloration morph") +
  theme(plot.title = element_text(hjust = 0.5))

