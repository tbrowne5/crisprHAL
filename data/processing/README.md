# Data processing for raw read files

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699

NCBI SRA Bioproject: PRJNA939699

* 1. Merge paired read files available at the SRA link
* 2. Separate sgRNAs from the merged reads file by their barcode using the perl script
* 3. Generate a file of 28 nucleotide inputs and 20 nucleotide sgRNA target sites with the python program prepare_28nt_sequences.py
* 4. In R, use ALDEx2 CLR and effect functions on the dataset with 10 selective and 10 non-selective conditions (details found at https://bioconductor.org/packages/release/bioc/html/ALDEx2.html)
* 5. Extend the 20 nucleotide sgRNA target sites by 8 nucleotides downstream to obtain the model input sequences


3. Python generation of the 28 nucleotide inputs to be merged in R
```
python prepare_28nt_sequences.py [Input_fasta_sequence]
```

4. R data processing to generate scores from raw data
```
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

```
