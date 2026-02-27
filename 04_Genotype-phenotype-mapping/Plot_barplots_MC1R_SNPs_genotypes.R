## Barplots - Genotypes at MC1R SNPs:
## Capture dataset analysis
## Last edited: 20-10-2023, Inês Miranda


#Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


####Genotypes gwasCOAT dataset####

dataCOAT <- read.table("./genotypes_MC1R_SNPs_gwasCOAT_dataset.txt", 
                   header = TRUE, sep = "\t", na.strings = "N" )
str(dataCOAT)

#Produce summary table w/ barplot information:

bardataCOAT <- table (dataCOAT$Genotypes, dataCOAT$Morph)
bardataCOAT

barplot (bardataCOAT, xlab = "Winter colour",
         ylab = "No. of specimens",
         ylim = c(0, 70),
         col = c("grey90", "bisque3", "tan4"), beside = TRUE,
         legend = row.names(bardataCOAT),
         cex.axis = 1.25,
         cex.names = 1.5,
         cex.lab = 1.5)


###Genotypes gwasLINE dataset####

dataLINE <- read.table("./genotypes_MC1R_SNPs_gwasLINE_dataset.txt", 
                   header = TRUE, sep = "\t", na.strings = "N" )
str(dataLINE)

#Produce summary table w/ barplot information:

bardataLINE <- table (dataLINE$Genotypes, dataLINE$Morph)
bardataLINE

barplot (bardataLINE, xlab = "Dorsoventral line",
         ylab = "No. of specimens",
         ylim = c(0, 70),
         col = c("grey90", "bisque3", "tan4"), beside = TRUE,
         legend = row.names(bardataLINE),
         cex.axis = 1.25,
         cex.names = 1.5,
         cex.lab = 1.5)


###Genotypes gwasMORPHS dataset####

dataMORPHS <- read.table("./genotypes_MC1R_SNPs_gwasMORPHS_dataset.txt", 
                   header = TRUE, sep = "\t", na.strings = "N" )
str(dataMORPHS)

#Produce summary table w/ barplot information:

bardataMORPHS <- table (dataMORPHS$Genotypes, dataMORPHS$Morph)
bardataMORPHS

barplot (bardataMORPHS, ylab = "Colour morph",
         xlab = "No. of specimens",
         col = c("grey90", "bisque3", "tan4"), beside = TRUE,
         xlim = c(0,120),
         horiz = TRUE,
         legend = row.names(bardataMORPHS),
         cex.axis = 1.25,
         cex.names = 1.5,
         cex.lab = 1.5)

# Note that MC1R is encoded in the 3'-5' strand of the reference genome and thus 
# genotypes are the complement of those shown in Supplementary Table 9.
