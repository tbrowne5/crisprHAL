library(ALDEx2)

data <- read.csv("dataset_counts.txt",header=True,Sep="\t")
data <- t(data)
data <- data[,c(3,5,8,10,12,14,16,18,20,1,4,6,9,11,13,15,17,19,21,2)]

conds <- c(rep("NS",10),rep("S",10))
data.clr <- aldex.clr(data,conds)
data.effect <- aldex.effect(data.clr)
data.output <- data.effect

sgRNAs <- read.csv("sgRNA_20_to_28nt_fasta_name.csv",header=TRUE,sep="\t")
output <- merge(sgRNAs,data.effect,by.x=2,by.y=0)
output <- output[,c(2,6)]
row.names(output) <- output[,1]
output[,1] <- NULL
output[,1] <- output[,1] / sd(output[,1]

write.csv(output,"sgRNA_inputs_and_scores.csv",quote=FALSE)
