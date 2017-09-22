#Filter data based on fdr values (PS_analysis directory)

#table = read.table("codeml_bind_data_qvals.csv", header=TRUE, sep=",")
#head(table)
#plot(density(log(table$lfdr)), xlim=c(-80,-5), ylim=c(0,0.01))
#ans <- table %>% filter(log(lfdr) < -5)
#ans
#unique(ans$CAZy_Family)
#ans %>% group_by(CAZy_Family)
#filtered <- ans %>% group_by(CAZy_Family) %>% summarise(n = n())
#filtered

# Write out the data filtered based on certain fdr 
#final_filtered <- filter(final_data, final_data$lfdr <  )
#final_filtered
#write.csv(final_filtered, file = "final_filtered_data.csv")

# Filter final data based on qvalues  (PS_analysis directory)

table = read.table("no_external_nodes_ps.csv", header=TRUE, sep=",")
table
#plot(density(log(table$qvalue)), ylim=c(0,0.001))
#table %>% filter(log(qvalue) < -1)
plot(density(log(table$qvalue)), xlim=c(-80,10), ylim=c(0,0.001))
ans1 <- table %>% filter(log(qvalue) < -20)
ans1
unique(ans1$CAZy_Family)
ans1 %>% group_by(CAZy_Family)
filtered <- ans1 %>% group_by(CAZy_Family) %>% summarise(n = n())
head(ans1)
names(ans1) <- c("index","Unnamed..0", "X.3", "X.2", "X.1", "X", "CAZy_Family","tdir", "lnL0","lnL1","LRT","p.value", "lfdr", "qvalue", "localfdr", "identifier")

# Remove columns with group and index
ans1$Unnamed..0 <- NULL
ans1$X.3 <- NULL
ans1$X.2 <- NULL
ans1$X.1 <- NULL
ans1$X <- NULL
ans1$index <- NULL
head(ans1)
write.csv(filtered, file = "/Users/eirinixemantilotou/Documents/PhD/PS/filtered_cazy_fams.csv")
write.csv(ans1, file = "/Users/eirinixemantilotou/Documents/PhD/PS/final_ps.csv")
