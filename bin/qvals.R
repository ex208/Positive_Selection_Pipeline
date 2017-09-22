# Run function for fdr calculation
# call the function for multiple files within the qvals_2 directory 
# The directory stores several csv files one for each CAZy family including daya for fdr

source("/Users/eirinixemantilotou/Documents/PhD/PhD_Project/r_code/qvalue_function.R")

# Make directory in which we will store all csv files for split fams
mainDir <- "./"
subDir <- "qvals"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

files <- list.files("./")
print(files)
for (i in 1:length(files)) {
  qval_calc(read.table(files[i], header=TRUE, sep=","))
}

#
table %>% filter(log(lfdr) < -1)
# Filter final data based on qvalues  
final_data = read.table("../bind_data_qvals.csv", header=TRUE, sep=",")
plot(density(log(final_data$qvalue)), ylim=c(0,1))
table %>% filter(log(lfdr) < -1)

final_filtered <- filter(final_data, final_data$qvalues < 0.01 )
final_filtered
write.csv(final_filtered, file = "final_filtered_data.csv")

