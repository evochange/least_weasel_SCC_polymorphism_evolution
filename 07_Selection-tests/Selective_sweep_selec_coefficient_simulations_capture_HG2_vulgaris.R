## Selective sweep - selection coefficient estimate:
## Capture dataset analysis
## Last edited: 31-05-2024, Inês Miranda


## Set working directory:
setwd ("/path/to/folder")
rm(list = ls())


## Import alpha values from SweeD bootstraps:
alpha_vals <- scan("alpha-means.capture_homoBrown_HG2_scaf337.MC1R_region.txt") 

  ## Check alpha values:
  str(alpha_vals)
  length(alpha_vals)

  ## Plot alpha distribution:
  hist(alpha_vals)
  
  ## Get stats for alpha:
  min (alpha_vals)
  max (alpha_vals)
  quantile (alpha_vals, c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99))


## Estimate selection coefficient values:
s_vals = c()

## Loop over alpha values
## Other variables are sampled from pre-determined distributions
for (i in 1:1000) {
  a <- alpha_vals[i]
  r <- runif (1, min = 0.0000000056, max = 0.0000000126)
  Ne <- runif (1, min = 349140, max = 386651) #using Ne estimates from North ancestry
  ln <- log(2*Ne)
  s=(r*ln)/a
  s_vals [i] <- s
}

  ## Get stats for s values:
  hist (s_vals)
  mean (s_vals)

  ## Get s distribution:
  quantile (s_vals, c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99))

  ## Output s values:
  ## write (s_vals, "Selection_Coefficient_HG2-B_estimates_1k_replicates.txt", n=1)

  # ----------------
  
  ## Read in s_vals from the original simulation, reported in the paper:
  s_vals <- scan("./Selection_Coefficient_HG2-B_estimates_1k_replicates.txt")
  
  # ----------------
  
## Find most likely s value (from density distribution):
alpha_s <- data.frame(alpha_vals, s_vals)
str (alpha_s)
alpha_s <- alpha_s[alpha_s$s_vals > 0 ,]

  ## Get y value with max density:
  which.max(density(alpha_s$s_vals)$y)
  #[1] 149

  ## Get corresponding y value
  density(alpha_s$s_vals)$x[149]
  #[1] 0.003302197

  ## Produce S density plot:
  d <- density(alpha_s$s_vals)
  d


## Plot S Values density plot:

  ## Import library:
  library (ggplot2)

  ## Draw plot:
  ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + 
    theme_classic() +
    geom_density(stat = "identity", colour = '#8b5a2b', fill = '#8b5a2b', alpha =0.1, size = 0.6) +
    geom_vline(xintercept = density(alpha_s$s_vals)$x[149], lty = 2 , col = "gray40", size = 1) +
    scale_y_continuous(ylab ("Density"), breaks =c(0,75,150), limits = c(0,150)) +
    xlab ("Selection Coefficient (s)")

