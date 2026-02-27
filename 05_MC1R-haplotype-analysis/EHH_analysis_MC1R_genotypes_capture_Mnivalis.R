## Extended Haplotype Homozygosity analysis:
## Capture dataset analysis
## Last edited: 20-12-2023, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


## Import library:
## install.packages("rehh")
library (rehh)


##### EHH ANALYSES WITH WITH PHENOTYPED DATASET (Sample Set 1) #####

## Import VCF with Sample set 1 (phenotyped M. nivalis):
gwasB_vcf2hh <- data2haplohh(hap_file = "./capture_phenotypedMnivalis_mc1r_flt1.variants.PSphased.ano.vcf",
                         vcf_reader = "data.table",
                         polarize_vcf = FALSE)
str (gwasB_vcf2hh)


## Computation of (allele-specific) EHH
## for the detection of selection within a single, homogeneous, population

gwasB_res <- calc_ehh (gwasB_vcf2hh, mrk = "337:165924",
                       discard_integration_at_border = FALSE,
                       limehh = 0.005)
str (gwasB_res)


## View results:
gwasB_res

## Plot results for EHH per allele
plot (gwasB_res, col = c("grey50", "tan4"),
      lwd = 2, cex.axis = 1.2, cex.lab = 1.2)
## Note that "ancestral" in the plot is the nivalis allele (equal to the reference genome)
## whereas "derived" in the plot is the vulgaris allele (alternative to the reference)

