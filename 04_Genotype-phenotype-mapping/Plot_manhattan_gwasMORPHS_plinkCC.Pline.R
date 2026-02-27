## Plot Case Control - Color Morph:
## Capture dataset analysis
## Last edited: 19-10-2023, Inês Miranda

library(qqman)

tiff("capture_gwasMORPHS_capture_Mnivalis.flt1.CC.BONF.Pline05.tif", width = 2500, height = 600, units = "px", res = 150)

data <- read.table("capture_gwasMORPHS_capture_Mnivalis.flt1.CC.BONF.valid.2plot", header = TRUE)

manhattan(data, chr = "CHR", bp = "BP", p = "P", snp = "SNP", logp = FALSE, cex = 1, cex.axis = 1, cex.axis = 1, cex.lab = 1, ylim = c(0, 60), col = c("gray60", "gray10"), chrlabs = c(1:231), ylab = "CC/-log10(P)", suggestiveline = FALSE, genomewideline = 6.641, xlab = "Scaffolds")

