import sys, re,datetime
from Bio.Seq import Seq

outname = sys.argv[1].split(".")
inputfile = open(sys.argv[1],"r")
outpath = sys.argv[2] # Better to allow for user to fully define output location
outfile = open(outpath,"w+")
print("\nPreparing the file: " + outpath)
outfile.write("sgRNAs\n")

forward = 0
reverse = 0
dups = 0
seq = ""
seqs = [] # List to store DNA seq of each contig

for line in inputfile:
    if line[0] != ">": 
        seq = seq + line.strip("\n").upper()
    if line[0] == ">" and len(seq) > 0: 
        seqs.append(seq) # Add contig sequence to list
        seq="" # Clear contig sequence to store next contig sequence

def find_sequences(dna_string):
    # Use a regex pattern to find all occurrences of 'GG'
    pattern = re.compile(r'(?=([ATCG]{21}GG[ATCG]{5}))') # Takes 21 upstream and 5 downstream nucleotides of GG

    # Find all matches
    matches = pattern.finditer(dna_string)

    # Initialize an empty list to store the results
    sequences = []

    # Iterate over the matches to extract the sequences
    for match in matches:
        # Extract the sequence and add it to the list
        sequence = dna_string[match.start(1):match.end(1)]
        sequences.append(sequence)

    return sequences

full_seqs = [] # As you loop through contigs, add sgRNA seqs to this list
for seq in seqs:          
    revstr = Seq(seq).reverse_complement() # BioPython reverse complement for potential runtime boost
    seq = seq + seq[0:27]
    revstr = str(revstr) + str(revstr)[0:27]
    full_seqs = full_seqs +find_sequences(seq) # add forward strand sgRNA
    full_seqs = full_seqs +find_sequences(revstr) # add reverse strand sgRNA

dedup_seqs = list(set(full_seqs)) # by using set, you will deduplicate the sequences in one shot instead of testing at every loop iter
for seq in dedup_seqs:
    outfile.write(seq+"\n")
outfile.close()
