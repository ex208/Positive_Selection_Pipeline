# Load eeptools
library(eeptools)

# Open data

files <- list.files(pattern = "*qvalue_lfdr.csv")
for (i in 1:length(files)) {
  #print(files[i])
  p = read.table(files[i], header=TRUE, sep=",")
  data(p)
  p <- p[order(p$qvalue) , ]
  print(p)
  write.csv(p, paste(filename = "/Users/eirinixemantilotou/Documents/PhD/PS/qvals/",files[i],"_ordered_qvalue_lfdr.csv", sep=""))
}