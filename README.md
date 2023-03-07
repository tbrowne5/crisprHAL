# crisprHAL: A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets

The CRISPR/Cas9 nuclease from Streptococcus pyogenes SpCas9 can be used with single guide RNAs (sgRNAs) as a sequence-specific antimicrobial agent and as a genome-engineering tool. However, current bacterial sgRNA activity models poorly predict SpCas9/sgRNA activity and are not generalizable, possibly because the underlying datasets used to train the models do not accurately measure SpCas9/sgRNA cleavage activity and cannot distinguish cleavage activity from toxicity. We solved this problem by using a two-plasmid positive selection system to generate high-quality biologically-relevant data that more accurately reports on SpCas9/sgRNA cleavage activity and that separates activity from toxicity. We developed a new machine transfer learning architecture (crisprHAL) that can be trained on existing datasets and that shows marked improvements in sgRNA activity prediction accuracy when transfer learning is used with small amounts of high-quality data. The crisprHAL model recapitulates known SpCas9/sgRNA-target DNA interactions and provides a pathway to a generalizable sgRNA bacterial activity prediction tool.

# Requirements:
```
python==3.8.8
numpy=1.19.2
numpy-base==1.19.2
biopython==1.78
h5py==2.10.0
hdf5==1.10.6
keras-preprocessing==1.1.2
pandas==1.2.2
scikit-learn==0.24.1
scipy==1.6.1
tensorflow==2.4.0
```

# Run Model Test:
```
python crisprHAL.py
```

# Process nucleotide sequence(s) into sgRNA target model inputs:
```
python process_inputs.py [Input Nucleotide File]
```

Input Nucleotide File: Single line nucleotide input with multiple sequences broken up by ">NAME" lines
Example input file shown below:
```
>Sequence1
TCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCC
>Sequence2
CCGACTCGTCCAACATCAATACAACCTATTAATTTCCCCTCGTCAAAAATAAGGTTATCAAGTGAGAAATCACCATGAGTGAC
>Sequence3
CCAGACTCCTGTGTAACATATGCAACCGTTCTAACCCGCTGGGTGAAGACTTTGACTACCGCAAAGAGTTTAGCAAGTTAGACTACTCCGCCCTGAAAGGGGATC
```

# Run Model:
```
python crisprHAL.py [Enzyme] [Input File csv] [Optional Compare]
```

Enzyme: "TevSpCas9" or "SpCas9"

Input: CSV file input name (format shown below)

Optional Compare: "Compare" if your CSV file contains activity scores for comparison (format shown below)

28 nucleotide sgRNA target site input required:
• 20 nucleotide target site, ie: CTCGATTGAGGGGCTGGGAA
• 3 nucleotide NGG PAM, ie: TGG
• 5 nucleotides downstream, ie: 
• Total: CTCGATTGAGGGGCTGGGAATGGGTGAT

Example Command and Input CSV File -- Prediction Only, No "Compare" Option:
```
python crisprHAL.py TevSpCas9 test_dataset.csv
```
```
sgRNA
ATGCATATCCCTCTTATTGCCGGTCGCG
GTCTTTATCAGCTAACCAGTCGGTATCC
CGATGGTCAATATCAGCCGTTGGCGCAG
TTTGCCTCATCAACACCTGAAGGCCTCA
CCGCTTTTCCTGCCTGACCTTGGGTGAA
CATAATAGTATTTCCGATAAGGGTCCCC
CGTAGCCGACAATACGGAAACGGTGAGT
CCGAGACGTTGATGCCAATATGGAAATT
CAGGAAACGGCTAACAGAACCGGACCAA
GTGGCAATCGTCGTTTTAACCGGCAAAC
```

Example Command and Input CSV File -- "Compare" option to compare to other activity scores:
```
python crisprHAL.py SpCas9 test_dataset.csv Compare
```
```
sgRNA,score
CTCGATTGAGGGGCTGGGAATGGGTGAT,8.21062788839
ATCTTTATCGGTCTTAGCAAAGGCTTTG,30.1092205446
CGGGCCAGACTGGCCGAGACGGGTCGTT,11.0586722001
CAGCATCATGCCGTGATCGCCGGGAAAG,0.67958058668
GGGGAAAGGACAAGCAGTGCGGGCAAAA,29.2338752707
GAGAGAATTTTGACCTGATCCGGCTCGC,2.28311187737
CGTGATGCAACTGTGTAATGCGGCTGAC,27.8665335936
GAACATCACCGCCTCACGTCCGGTTTTG,26.3480104838
TCGATTGAGGGGCTGGGAATGGGTGATC,41.2590972746
CCGTGTAAGGGAGATTACACAGGCTAAG,4.25926295656
```
