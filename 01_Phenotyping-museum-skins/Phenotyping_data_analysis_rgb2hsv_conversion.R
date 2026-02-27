## Phenotyping dataset - museum skins
## RGB to HSV data conversion
## Last edited: 26-05-2022, Inês Miranda


## Set working directory:
setwd("/pth/to/folder")
rm(list=ls())


## Read input:
m <- read.table("./20220522_NRM_museum_skins_phenotyping_RGB_raw_data.txt", 
                header=T, sep = "\t", na.strings = "n/a")
str (m)
table (m$MORPH, m$SEASON)

## Remove rows (individuals) with missing data:
rgb_indC <- m[complete.cases(m), ]
complete.cases (rgb_indC)
str(rgb_indC)

## Info on number of specimens remaining per category:
table(rgb_indC$SEASON,rgb_indC$MORPH)

## Transpose rows to columns:
head (rgb_indC)

rgb_indCt <- as.data.frame(t(rgb_indC))
head (rgb_indCt)
str(rgb_indCt)


## Create independent matrices for each sampling point:

LD_rgb <- as.matrix((rgb_indCt[c("LD.Red", "LD.Green", "LD.Blue"),]))
LD_rgb <- matrix (as.numeric(LD_rgb), ncol = ncol (LD_rgb))
row.names(LD_rgb) <- c("red", "green", "blue")
str (LD_rgb)

MD_rgb <- as.matrix((rgb_indCt[c("MD.Red", "MD.Green", "MD.Blue"),]))
MD_rgb <- matrix (as.numeric(MD_rgb), ncol = ncol (MD_rgb))
row.names(MD_rgb) <- c("red", "green", "blue")
str (MD_rgb)

UD_rgb <- as.matrix((rgb_indCt[c("UD.Red", "UD.Green", "UD.Blue"),]))
UD_rgb <- matrix (as.numeric(UD_rgb), ncol = ncol (UD_rgb))
row.names(UD_rgb) <- c("red", "green", "blue")
str (UD_rgb)

H_rgb <- as.matrix((rgb_indCt[c("H.Red", "H.Green", "H.Blue"),]))
H_rgb <- matrix (as.numeric(H_rgb), ncol = ncol (H_rgb))
row.names(H_rgb) <- c("red", "green", "blue")
str (H_rgb)

MV_rgb <- as.matrix((rgb_indCt[c("MV.Red", "MV.Green", "MV.Blue"),]))
MV_rgb <- matrix (as.numeric(MV_rgb), ncol = ncol (MV_rgb))
row.names(MV_rgb) <- c("red", "green", "blue")
str (MV_rgb)

UV_rgb <- as.matrix((rgb_indCt[c("UV.Red", "UV.Green", "UV.Blue"),]))
UV_rgb <- matrix (as.numeric(UV_rgb), ncol = ncol (UV_rgb))
row.names(UV_rgb) <- c("red", "green", "blue")
str (UV_rgb)

Hv_rgb <- as.matrix((rgb_indCt[c("H_v.Red", "H_v.Green", "H_v.Blue"),]))
Hv_rgb <- matrix (as.numeric(Hv_rgb), ncol = ncol (Hv_rgb))
row.names(Hv_rgb) <- c("red", "green", "blue")
str (Hv_rgb)


## Transform RGB values to HSV per matrix:

LD_hsv <- rgb2hsv(LD_rgb, maxColorValue = 65535)
row.names (LD_hsv) <- c("LD_h", "LD_s", "LD_v")
str (LD_hsv)

MD_hsv <- rgb2hsv(MD_rgb, maxColorValue = 65535)
row.names (MD_hsv) <- c("MD_h", "MD_s", "MD_v")
str (MD_hsv)

UD_hsv <- rgb2hsv(UD_rgb, maxColorValue = 65535)
row.names (UD_hsv) <- c("UD_h", "UD_s", "UD_v")
str (UD_hsv)

H_hsv <- rgb2hsv(H_rgb, maxColorValue = 65535)
row.names (H_hsv) <- c("H_h", "H_s", "H_v")
str (H_hsv)


## Binding all matrices:
total_hsv <- rbind (LD_hsv, MD_hsv, UD_hsv, H_hsv)

## Recover specimens data:
indv = rgb_indC[,c(1:3)]
indv$Group <- paste (indv$MORPH, indv$SEASON)
str (indv)

## Binding all info together:
total_hsv_t <- t(total_hsv)
totaldataset_hsv <- cbind(indv, total_hsv_t)
str(totaldataset_hsv)

## Write output:
write.table(totaldataset_hsv, file = "HSV_transformation_final-dataset_dorsal_data_27-05-2022.txt",
            quote = F, sep = "\t", row.names = F)
