# Make directory in which we will store all csv files for split fams
mainDir <- "./"
subDir <- "fdr"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)


source("/Users/eirinixemantilotou/Documents/PhD/PhD_Project/PS_analysis/r_code/fdr_function.R")

# Run function for fdr calculation
# call the function for multiple files within the split_fams directory 
# The directory stores several csv files. One csv for each CAZy family

files <- list.files("./")
for (i in 1:length(files)) {
  fdr_calc(read.table(files[i], header=TRUE, sep=","))
}

# After calculating the fdr for each individual family we combine those csvs into one again and name it as 
# bind_data_fdr.csv using the **bind.R** code

#Filter data based on fdr values (PS_analysis directory)
table = read.table("codeml_bind_data_qvals.csv", header=TRUE, sep=",")
head(table)
plot(density(log(table$lfdr)), xlim=c(-80,-5), ylim=c(0,0.01))
ans <- table %>% filter(log(lfdr) < -5)
ans
#unique(ans)
unique(ans$CAZy_Family)
ans %>% group_by(CAZy_Family)
filtered <- ans %>% group_by(CAZy_Family) %>% summarise(n = n())
filtered

# Write out the data filtered based on certain fdr 
#final_filtered <- filter(final_data, final_data$lfdr <  )
#final_filtered
#write.csv(final_filtered, file = "final_filtered_data.csv")