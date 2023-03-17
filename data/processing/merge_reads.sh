################################################################################
# Insert the path to your usearch install and downloaded reads.
usearch='/Volumes/data/bin/usearch'
READS='/Volumes/Downloads/crisprHAL/readprocessing/reads/'
################################################################################

# make output directory for usearch files
mkdir -p usearch_output

for f in ./split_reads/*R1*.fastq; do

#get basename
B=`basename $f`
NAME=`echo $B | cut -d "." -f1`
SAMPLE=`echo $B | cut -d "_" -f1`

R='R2'
# replace R1 with R2 to specify reverse read file
REV=`echo ${B/R1/$R}`
REVNAME=`echo ${NAME/R1/$R}`

$usearch -report ./usearch_output/$SAMPLE.report -fastq_mergepairs ./reads/${SAMPLE}_R1.fastq -reverse ./reads/${SAMPLE}_R2.fastq -fastqout ./usearch_output/$SAMPLE.fastq -fastaout ./usearch_output/$SAMPLE.fasta

echo $SAMPLE
echo $SAMPLE.done

done