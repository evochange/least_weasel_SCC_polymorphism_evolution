## Plot TreeMix:
## Capture dataset analysis
## Last edited: 06-10-2025, Inês Miranda

## Set working directory
setwd ("/path/to/file")
rm(list = ls())

## Import plotting functions file from TreeMix:
## File available at: https://bitbucket.org/nygcresearch/treemix/src/master/src/
source("./plotting_funcs.R")

## Plot tree for migration edge value =  0:
plot_tree("treemix_europe3ANC_1pop_outgp_2k_frag_regions.mac2.25k.mig0", cex=0.8)

## Plot residuals for the  tree:
plot_resid("treemix_europe3ANC_1pop_outgp_2k_frag_regions.mac2.25k.mig0", "./pops-order.txt", cex=0.7)
