# Data processing for raw read files

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699

NCBI SRA Bioproject: PRJNA939699

* 1. Merge paired read files available at the SRA link
* 2. Separate sgRNAs from the merged reads file by their barcode using the perl script
* 3. Load the resulting file into R
* 4. Run the R script with ALDEx2 to generate the standard deviation scaled difference-between value scores
