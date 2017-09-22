# Code for generating q values based on given p values, plots for q values and p values
library(dplyr) 
library(qvalue)

# Function to calculate qvalues based on prefixed fdr and pvalues
qval_calc <- function(x) {
  #hist(x$V7, breaks = 20, main = paste("Distribution of p-values"), xlab="Value")
  qobj <- qvalue(x$p.value, pi0.meth="bootstrap", fdr.level=0.05)
  qvalues <- qobj$qvalues
  localFDR <- qobj$lfdr
  x["qvalue"] <- qvalues
  x["localfdr"] <- localFDR 
  write.csv(x, paste(file = "/Users/eirinixemantilotou/Documents/PhD/PS/qvals/",files[i],"_qval.csv", sep=""))
}

