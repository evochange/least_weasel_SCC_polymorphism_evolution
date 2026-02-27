## MC1R haplotype heatmap (phased data):
## Capture dataset analysis
## Last edited: 02-02-2024, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


## Load libraries:
library(tidyverse)
library (dplyr)
library(ggplot2)


#### IMPORT DATA ####

## Read input files
table<-read.table("./capture_europeMnivalis_mc1r_flt1.variants.maf10.snps40kb.PSphased.haplo")
str(table)

## Read the header from the vcf (remove everything but CHR, POS, and genotype columns)
header <- read.table("./capture_europe_mc1r_flt1.haplotypes.header.txt", sep = "\t")
str (header)

## Add header to haplotypes table:
colnames(table)<-header[,1]
str(table)

## Write new table to file:
str(table)
write.table (table,"./capture_europeMnivalis_mc1r_flt1.variants.maf10.snps40kb.PSphased.header.haplo", sep ="\t")


#### PLOTING ####

## Read the list of individuals ordered as you want them to be plotted:
indv<-read.table("./capture_europe_mc1r_flt1.samples_order.txt")
str (indv)


## Convert genotype table into a longer format:
table %>% 
  pivot_longer(cols=3:486,names_to="individuals",values_to="genotypes") -> table_genos
str(table_genos)

table_genos$POS<-as.character(table_genos$POS)
table_genos$genotypes<-as.character(table_genos$genotypes)
str(table_genos)


## Choose color ramp for plotting: 
colors <- c("gray95", "grey20")


## Plot haplotype heatmap:
p1<- table_genos %>%
  filter(CHROM=="337") %>%
  
  ggplot() +
  geom_tile(aes(x=POS,y=factor(individuals,level=indv[,1]),fill=genotypes))+
  scale_fill_manual(values=colors)+
  xlab("Individuals") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
p1

## Save heatmap:
ggsave("./capture_europeMnivalis_mc1r_flt1.variants.maf10.snps40kb.PSphased.haplotypes.tiff", units="in", width=40, height=20, dpi=300, limitsize=FALSE, compression = 'lzw')

