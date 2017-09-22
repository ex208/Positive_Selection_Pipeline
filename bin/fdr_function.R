# Code for generating local fdr values based on given p values, plots for fdr and p values
library(dplyr)
library(ggplot2)
library(fdrtool)


# Function to calculate fdr values based on pvalues
fdr_calc <- function(x) {
  #x <- x %>% filter(is.na(p.value ))
  pval <- (x$p.value)
  fdr = fdrtool(pval, statistic="pvalue")
  x["lfdr"]<- fdr$lfdr
  
  #plot(density(log(p$fdr)), ylim=c(0,100000))
  #p %>% filter(log(fdr) < -1)
  write.csv(x, paste(file = "/Users/eirinixemantilotou/Documents/PhD/PS/fdr/",files[i],"_fdr.csv", sep=""))
}


