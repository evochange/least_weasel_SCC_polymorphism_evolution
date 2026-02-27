## MC1R haplogroup counts - barplots:
## Capture dataset analysis
## Last edited: 05-02-2024, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


####Barplots per HAPLOGROUP####

data <- read.table("./capture_europeMnivalis_MC1R_haplogroups_ancestry_information.txt", 
                   header = TRUE, sep = "\t")
str(data)

#Produce summary table w/ barplot information:

bardata <- table (data$Haplotype_Group, data$Admixture_Ancestry)
bardata

barplot (bardata, xlab = "Population Ancestry",
         ylab = "No. of haplotypes",
         ylim = c(0, 160),
         col = c("#8b5a2b", "#9F9F9F", "#cdb79e", "#e5e5e5", "#ffffff"),
         #density = c(30, 30, 20, 20, 20),
         #angle = c(45, 45, 90, 90, 0),
         beside = TRUE,
         legend = row.names(bardata),
         cex.axis = 1.25,
         cex.names = 1.5,
         cex.lab = 1.5)
#Note that "B" stand for vulgaris allele and "W" stands for nivalis allele
