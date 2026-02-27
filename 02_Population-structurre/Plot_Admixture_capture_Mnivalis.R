## Plot ADMXITURE:
## Capture dataset analysis
## Last edited: 06-10-2025, Inês Miranda

## Original author: Joana Meier, September 2019
## Original code from: https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r


#Set working directory:
setwd ("/path/to/folder")
rm(list=ls())


#### PLOT CROSS VALIDATION ERROR ####

## Read input data:
cv_data <- read.table ("./CV_error_capture_europe_genomic_2k_snps_regions_mac2.dist25k.txt", sep = " ")
head (cv_data)

## Add column names:
colnames(cv_data) <- c("K_value", "CV")
str (cv_data)

## Plot results:
library (ggplot2)

ggplot (cv_data, aes (x=K_value, y=CV)) +
  theme_bw() +
  geom_line(colour = "black", linewidth = 1.2) + 
  xlab ("No. ancestral culster (K)") +
  ylab ("Cross-validation error") +
  theme(axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size=14),
      axis.text.y = element_text(size=14))
      

#### PLOT ADMIXTURE RESULTS ####

## Assign the file prefix name:
prefix <- "capture_europe_genomic_2k_snps_regions_mac2.dist25k"

## Get individual names in the correct order:
labels<-read.table("capture_europe_3lineages_individual_country_info.txt")

## Name the columns:
names(labels)<-c("ind","pop")
str(labels)

## Define population order:
populations <- "FIN,SWE,SWE2,DNK,LUX,SWI,GER,POL,AUT,GER2,POL2,AUT2,CZE"

## Define individual order:
individuals <- readLines("./capture_europe_3lineages_reordering_individuals.csv")
individuals

## Add a column with population indices to order the barplots:
## Use the order of populations provided
labels$n<-factor(labels$pop,levels=unlist(strsplit(populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))
str(labels)

## Add a column with individual indices to order the barplots:
## Use the order of individuals provided
labels$nind<-factor(labels$ind,levels=unlist(strsplit(individuals,",")))
levels(labels$nind)<-c(1:length(levels(labels$nind)))
labels$nind<-as.integer(as.character(labels$nind))
str(labels)


## Read in the different Admixture output files:
minK=1  #note: a fake assignment for K = 1 was created to include the phenotypic assignment of each sample in the plot
maxK=3
tbl<-lapply(minK:maxK, function(x) read.table(paste0("./admixture_output_files/",prefix,".",x,".Q")))

## Prepare spaces to separate the populations:
rep<-as.vector(table(labels$n))
rep
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
spaces<-spaces[-length(spaces)]


## Plot the cluster assignments as a single bar for each individual, with each K as a separate row

par(mfrow=c(maxK,1),mar=c(0,1,0,0),oma=c(2,1,1,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)

#note: this first line will show phenotypic assigment (nivalis, vulgaris, uncertain) for each sample
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$nind),])), col=c("white","tan4","grey"),xaxt="n", border=NA,ylab=paste0("K=",minK),yaxt="n",space=spaces)

bp <- barplot(t(as.matrix(tbl[[2]][order(labels$nind),])), col=c("#DDAA33","#004488"),xaxt="n", border=NA,ylab=paste0("K=",minK+1),yaxt="n",space=spaces)

barplot(t(as.matrix(tbl[[3]][order(labels$nind),])), col=c("#DDAA33","#BB5566", "#004488"),xaxt="n", border=NA,ylab=paste0("K=",minK+2),yaxt="n",space=spaces)

axis(1,at=c(which(spaces==0.5),bp[length(bp)])-diff(c(1,which(spaces==0.5),bp[length(bp)]))/2,
     labels=unlist(strsplit(populations,",")))

