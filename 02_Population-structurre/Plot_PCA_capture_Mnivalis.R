## Plot Principal Component Analysis from PLINK:
## Capture dataset analysis
## Last edited: 06-10-2025, Inês Miranda

## load libraries:
require(optparse)
require(ggplot2)
require(ggrepel)

## Set working directory:
setwd("/path/to/folder")
rm(list=ls())

## Read PLINK output with eigenvectors:
## Output edited to include information on samples and their morph and assigned genetic lineage
pca_eigenvec <- read.table("./capture_europe_genomic_2k_snps_regions_mac2.dist25k.eigenvec.txt", sep = "\t", header=T)
str (pca_eigenvec)

## Read PLINK output with eigenvalues:
pca_eigenval <- scan ("./capture_europe_genomic_2k_snps_regions_mac2.dist25k.eigenval")
pca_eigenval

## Calculate percentage of explained variance in each PC:
pev <- data.frame (PC = 1:20, pev = pca_eigenval/sum(pca_eigenval)*100)
pev

## Create plot with percentages of explained variance:
ggplot(pev, aes (PC,pev)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  ylab ("Percentage explained variance")

## Calculate cumulative sum of explained variance:
cumsum(pev$pev)

## Plot PCA graph:
a <- ggplot(pca_eigenvec, aes(x=(PC1*-1), y=PC2, col=Population, shape=Morph, fill=Morph)) +  #multiply PC1 by -1 to have PCA groups in the same order as Admixture
  theme_bw() +
  geom_point(size = 4, stroke = 1)
a <- a + scale_shape_manual(name = "Morph", values = c(21,23,22)) 
a <- a + scale_color_manual(values = c("#DDAA33","#BB5566","#004488","#BBBBBB"))
a <- a + scale_fill_manual(values = c("#ffffff", "#bebebe", "#8b5a2b"))
a <- a + xlab(paste0("PC1 (", round(pev$pev[1],2),"%)")) +
  ylab(paste0("PC2 (", round(pev$pev[2],2),"%)"))+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 16))
a
