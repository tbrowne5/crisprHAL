# crisprHAL: A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets

The CRISPR/Cas9 nuclease from Streptococcus pyogenes SpCas9 can be used with single guide RNAs (sgRNAs) as a sequence-specific antimicrobial agent and as a genome-engineering tool. However, current bacterial sgRNA activity models poorly predict SpCas9/sgRNA activity and are not generalizable, possibly because the underlying datasets used to train the models do not accurately measure SpCas9/sgRNA cleavage activity and cannot distinguish cleavage activity from toxicity. We solved this problem by using a two-plasmid positive selection system to generate high-quality biologically-relevant data that more accurately reports on SpCas9/sgRNA cleavage activity and that separates activity from toxicity. We developed a new machine transfer learning architecture (crisprHAL) that can be trained on existing datasets and that shows marked improvements in sgRNA activity prediction accuracy when transfer learning is used with small amounts of high-quality data. The crisprHAL model recapitulates known SpCas9/sgRNA-target DNA interactions and provides a pathway to a generalizable sgRNA bacterial activity prediction tool.

# Sections of this guide:

Setting up and running the model to predict sgRNA activities:
* 0: Model requirements
* 1: Running the model test
* 2: Processing nucleotide sequences into model inputs
* 3: Predicting with the model

Additional information and methods: 
* 4: Preparing your own model input files & comparing predictions
* 5: Validating the trained models
* 6: Citations

# 0: Requirements:
```
python==3.8.8
numpy==1.19.2
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

# 1: Run model test:
```
python crisprHAL.py
```
Test our TevSpCas9 model with an example SpCas9 dataset of 7821 unique sgRNA target sites from Guo et al. 2018.



# 2: Process nucleotide sequence(s) into sgRNA target model inputs:

This will take an input nucleotide fasta file and identify potential sgRNA sequences for evaluation. The output
will be a .csv file containing the predicted guides. This can be used as input for the prediction step

```
# python process_inputs.py [Input Nucleotide File]
python process_inputs.py phiX174.fna

#output phiX174_output.csv
```

* **Input Nucleotide File**: Single line nucleotide input with multiple sequences broken up by ">NAME" lines
* **Output**: 28 nucleotide sequences in a CSV file appropriate as an input to the model

Composition of the 28 nucleotide inputs:
* 20 nucleotide target site, ie: CTCGATTGAGGGGCTGGGAA
* 3 nucleotide NGG PAM, ie: TGG
* 5 nucleotides downstream, ie: GTGAT
* Total: CTCGATTGAGGGGCTGGGAATGGGTGAT

Example input file shown below:
```
>Sequence1
TCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCC
>Sequence2
CCGACTCGTCCAACATCAATACAACCTATTAATTTCCCCTCGTCAAAAATAAGGTTATCAAGTGAGAAATCACCATGAGTGAC
>Sequence3
CCAGACTCCTGTGTAACATATGCAACCGTTCTAACCCGCTGGGTGAAGACTTTGACTACCGCAAAGAGTTTAGCAAGTTAGACTACTCCGCCCTGAAAGGGGATC
```


# 3: Predict with model: 

This will take the file of the predicted sgRNA sequences from above and assign a score. Higher scores are better!
The output a .csv file named OUTPUT_[inputfile] and contains the sgRNA and the score.
```
# python crisprHAL.py [Enzyme] [Input file csv] [Optional compare]
python crisprHAL.py TevSpCas9  phiX174_output.csv
# output
```

* **Enzyme**: "TevSpCas9" or "SpCas9"
* **Input**: CSV file input name; created in section 2 or matching the required format (Format: Section 4)
* **Optional Compare**: "Compare" if your CSV file contains scores for comparison to model predictions (Format: Section 4)
* **Output**: Tab-deliminated (TSV) file containing the 28 nucleotide sequence and predicted Cas9 activity value

Example command with prediction only, no "compare" option:
```
python crisprHAL.py TevSpCas9 test_dataset.csv
```

Example command with prediction and "compare" option for prediction comparison:
```
python crisprHAL.py SpCas9 test_dataset.csv compare
```


# 4: Preparing your own input CSV Files:

Input CSV file for prediction only, no "Compare" option:
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

Input CSV file for prediction and comparison of the predictions to scores in column 2:
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

# 5: Validate the training of the model

Perform 5-fold cross validation with the TevSpCas9 dataset transfer learning from the eSpCas9 base model:
```
python crisprHAL.py train TevSpCas9
```

Perform 5-fold cross validation with the SpCas9 dataset transfer learning from the eSpCas9 base model:
```
python crisprHAL.py train SpCas9
```

Perform 80:20 train-test split with the eSpCas9 dataset (Guo et al. 2018) used as the base model\*:
```
python crisprHAL.py train eSpCas9
```
\*An 80% training & 20% test split was used for base model generation, and therefore has been included in place of the 5-fold cross validation tests used for the TevSpCas9 and SpCas9 enzyme transfer learning-based models.

# 6: Citations -- TO DO
