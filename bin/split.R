#This code has as purpose to split the initial file with all cazy results in individual cazy families

# Make directory in which we will store all csv files for split fams
mainDir <- "./"
subDir <- "split_fams"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

table = read.table("codeml_bind_data.csv", header=TRUE, sep=",")
head(table)
split_list <- split(table,table$CAZy_Family)
lapply(names(split_list), function(x){write.csv(split_list[[x]], file = paste( x,".csv", sep = ""))})

# Move all csv files in a directory called split_fams