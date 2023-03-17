# Data processing for raw read files

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699

NCBI SRA Bioproject: PRJNA939699

* 1. Merge paired read files available at the SRA link
* 2. Separate sgRNAs from the merged reads file by their barcode using the perl script
* 3. Load the resulting file into R
* 4. Use ALDEx2 CLR and effect functions on the dataset with 10 selective and 10 non-selective conditions (details found at https://bioconductor.org/packages/release/bioc/html/ALDEx2.html)
* 5. Extend the 20 nucleotide sgRNA target sites by 8 nucleotides downstream to obtain the model input sequences

R data processing to generate scores from raw data
```
library(ALDEx2)

data <- read.csv("dataset_counts.txt",header=True,Sep="\t")
data <- t(data)
data <- data[,c(3,5,8,10,12,14,16,18,20,1,4,6,9,11,13,15,17,19,21,2)]

conds <- c(rep("NS",10),rep("S",10))
data.clr <- aldex.clr(data,conds)
data.effect <- aldex.effect(data.clr)
data.output <- data.effect
data.output[,c(1,2,3,5,6,7)] <- NULL
data.output$diff.btw <- data.output$diff.btw / sd(data.output$diff.btw)
```
