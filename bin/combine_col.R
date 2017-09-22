# Code for combining names from two columns and excluding other columns
library(dplyr)
p = read.table("codeml_results_split.csv", header=FALSE, sep=",")
names(p) <- c("index","CAZy_Family","lnL0","group","tdir","lnL1","LRT","p.value")
head(p)
new_p = transform(p, CAZy_Family=paste(CAZy_Family, group, sep="_"))
# Remove columns with group and index
new_p$group <- NULL
new_p$index <- NULL
head(new_p)
new_p <- new_p[c("CAZy_Family", "tdir", "lnL0", "lnL1", "LRT", "p.value")]
head(new_p)
write.csv(new_p, file = "codeml_results_split_COMB.csv")
