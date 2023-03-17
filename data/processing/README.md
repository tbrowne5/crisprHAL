# Data processing for raw read files

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699

NCBI SRA Bioproject: PRJNA939699

Order of operations for file generation:
* 1. Merge the paired read files available at the SRA link
* 2. Separate sgRNAs from the merged reads file by the barcodes found in the file *primers.txt* (example perl script)
* 3. Generate a file of 28 nucleotide inputs and 20 nucleotide sgRNA target sites with the command: **python prepare_28nt_sequences.py [input_fasta_sequence]**
* 4. In R, use ALDEx2 CLR and effect functions on the dataset with the R script *get_inputs_and_scores.R* (ALDEx2 information available at https://bioconductor.org/packages/release/bioc/html/ALDEx2.html)
