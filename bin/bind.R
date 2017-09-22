# Make codeml_results csv same with codeml_results_split_COMB
results = read.table("codeml_results.csv", header=FALSE, sep=",")
names(results) <- c("index","CAZy_Family","lnL0","tdir","lnL1","LRT","p.value")
results$index <- NULL
new_results <- results[c("CAZy_Family", "tdir", "lnL0", "lnL1", "LRT", "p.value")]
head(new_results)
write.csv(new_results, file = "codeml_results_COMB.csv")


# 1. we combine the 2 csvs from the two separet analysis 

first_analysis <- read.csv("codeml_results.csv",header=T,sep=",")
second_analysis <- read.csv("codeml_results_split_COMB.csv",header=T,sep=",")
final_data <- rbind(first_analysis, second_analysis)
final_data
write.csv(final_data, file = "codeml_bind_data.csv")

# 2. Combine all fdr csv's from PS_analysis/fdr
files  <- list.files(pattern = '\\.csv')
tables <- lapply(files, read.csv, header = TRUE)
combined.df <- do.call(rbind , tables)
write.csv(combined.df, file = "codeml_bind_data_fdr.csv")

# 3. Combine all qvals csv's from PS/anaysis/qvals
files  <- list.files(pattern = '\\.csv')
tables <- lapply(files, read.csv, header = TRUE)
combined.df <- do.call(rbind , tables)
write.csv(combined.df, file = "codeml_bind_data_qvals.csv")