# Make directory in which we will store all csv files for split fams
mainDir <- "./"
subDir <- "qvals"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

# call the function for qvals calculation

source("/Users/eirinixemantilotou/Documents/PhD/PhD_Project/PS_analysis/r_code/qvalue_function.R")
# Run function for qvals calculation
# call the function for multiple files within the fdr directory 
# The directory stores several csv files. One for each CAZy family
files <- list.files("./")
print(files)
for (i in 1:length(files)) {
  qval_calc(read.table(files[i], header=TRUE, sep=","))
}


#Combine all data in one csv form qvals_2 directory using the bind.R code

# Filter final data based on qvalues  (PS_analysis directory) 
final_data = read.table("codeml_bind_data_qvals.csv", header=TRUE, sep=",")
#plot(density(log(final_data$qvalue)), ylim=c(0,0.001))
#table %>% filter(log(qvalue) < -1)
plot(density(log(table$qvalue)), xlim=c(-80,10), ylim=c(0,0.001))
ans <- table %>% filter(log(qvalue) < )
unique(ans$CAZy_Family)
ans %>% group_by(CAZy_Family)
filtered <- ans %>% group_by(CAZy_Family) %>% summarise(n = n())
filtered



# Write out the data filtered based on certain qvalue
#final_filtered <- filter(final_data, final_data$qvalues <  )
#final_filtered
#write.csv(final_filtered, file = "final_filtered_data.csv")

