## Plot Case Control - Zoom-ins:
## Capture dataset analysis
## Last edited: 24-10-2023, Inês Miranda


#Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


# Load libraries:
require (ggplot2)


#### A. CASE-CONTROL LINE ####

  ## Import case-control results:
  cc.line.scaf337 <- read.table("./capture_gwasLINE_capture_Mnivalis.flt1.CC.BONF.scaf337.valid.2plot", header = T)
  str (cc.line.scaf337)

  summary(cc.line.scaf337$P)


  ## Get outliers CC:
  cc.line.highP <- cc.line.scaf337 [cc.line.scaf337$P >= 30,]
  cc.line.mediumP <- cc.line.scaf337 [cc.line.scaf337$P < 30 & cc.line.scaf337$P >= 6.612,]
  cc.line.lowP <- cc.line.scaf337 [cc.line.scaf337$P < 6.612,]


  ## Plot results:

  ggplot (data = cc.line.lowP, mapping = aes (x=BP, y=P)) + 
    theme_classic() +
    geom_point (data=cc.line.lowP, colour='gray70', shape = 16, size = 2) +
    geom_point (data=cc.line.highP, colour='darkorange', shape = 23, size = 3, fill='darkorange')+
    geom_point (data=cc.line.mediumP, colour='gray20', shape = 16, size = 2)+
    scale_x_continuous(xlab("Scaffold 337"), limits= c (0, 500000), breaks = c(0, 100000, 200000,300000, 400000, 500000)) + 
    scale_y_continuous(ylab("CC -log10 (P)"), limits = c(0,40), breaks = c(0, 20, 40, 60, 80, 100)) +
    geom_segment(aes(x=0, xend=500000, y=6.612, yend=6.612), colour="#cc3311",linewidth=1, lty= "dashed") +
    theme(
      axis.title.x = element_text(size=13.5),
      axis.title.y = element_text(size=13.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    ) 


  
#### B. CASE-CONTROL COAT ####

  ## Import case-control results:
  cc.coat.scaf337 <- read.table("./capture_gwasCOAT_capture_Mnivalis.flt1.CC.BONF.scaf337.valid.2plot", header = T)
  str (cc.coat.scaf337)

  summary(cc.coat.scaf337$P)


  ## Get outliers CC:
  cc.coat.highP <- cc.coat.scaf337 [cc.coat.scaf337$P >= 35,]
  cc.coat.mediumP <- cc.coat.scaf337 [cc.coat.scaf337$P < 35 & cc.coat.scaf337$P >= 6.495,]
  cc.coat.lowP <- cc.coat.scaf337 [cc.coat.scaf337$P < 6.495,]


  #Plot results:

  ggplot (data = cc.coat.lowP, mapping = aes (x=BP, y=P)) + 
    theme_classic() +
    geom_point (data=cc.coat.lowP, colour='gray70', shape = 16, size = 2) +
    geom_point (data=cc.coat.highP, colour='darkorange', shape = 23, size = 3, fill='darkorange')+
    geom_point (data=cc.coat.mediumP, colour='gray20', shape = 16, size = 2)+
    scale_x_continuous(xlab("Scaffold 337"), limits= c (0, 500000), breaks = c(0, 100000, 200000,300000, 400000, 500000)) + 
    scale_y_continuous(ylab("CC -log10 (P)"), limits = c(0,40), breaks = c(0, 20, 40, 60, 80, 100)) +
    geom_segment(aes(x=0, xend=500000, y=6.495, yend=6.495), colour="#cc3311",linewidth=1, lty= "dashed") +
    theme(
      axis.title.x = element_text(size=13.5),
      axis.title.y = element_text(size=13.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    ) 


  
#### C. CASE-CONTROL MORPHSL ####

  ## Import case-control results:
  cc.morphs.scaf337 <- read.table("./capture_gwasMORPHS_capture_Mnivalis.flt1.CC.BONF.scaf337.valid.2plot", header = T)
  str (cc.morphs.scaf337)

  summary(cc.morphs.scaf337$P)


  ## Get outliers Fst:
  cc.morphs.highP <- cc.morphs.scaf337 [cc.morphs.scaf337$P >= 55,]
  cc.morphs.mediumP <- cc.morphs.scaf337 [cc.morphs.scaf337$P < 55 & cc.morphs.scaf337$P >= 6.641,]
  cc.morphs.lowP <- cc.morphs.scaf337 [cc.morphs.scaf337$P < 6.641,]


  ## Plot results:

  ggplot (data = cc.morphs.lowP, mapping = aes (x=BP, y=P)) + 
    theme_classic() +
    geom_point (data=cc.morphs.lowP, colour='gray70', shape = 16, size = 2) +
    geom_point (data=cc.morphs.highP, colour='darkorange', shape = 23, size = 3, fill='darkorange')+
    geom_point (data=cc.morphs.mediumP, colour='gray20', shape = 16, size = 2)+
    scale_x_continuous(xlab("Scaffold 337"), limits= c (0, 500000), breaks = c(0, 250000, 500000), expand = c(0.02, 0.02)) + 
    scale_y_continuous(ylab("CC -log10 (P)"), limits = c(0,60), breaks = c(0, 15, 30, 45, 60), expand = c(0.02, 0.02)) +
    geom_segment(aes(x=0, xend=500000, y=6.641, yend=6.641), colour="#cc3311",linewidth=1, lty= "dashed") +
    theme(
      axis.title.x = element_text(size=13.5),
      axis.title.y = element_text(size=13.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    ) 

