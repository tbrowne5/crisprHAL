## crisprHAL: A suite of bacterial sgRNA/Cas9 activity prediction models
Prediction models for bacterial SpCas9, TevSpCas9, eSpCas9, and SaCas9 nuclease activity prediction.

## The crisprHAL series papers, paper repositories, and web tool:
#### [The up-to-date crisprHAL repository](https://github.com/tbrowne5/crisprHAL)
#### [The online crisprHAL prediction tool](https://crisprhal.streamlit.app/)
#### [The original crisprHAL paper repository](https://github.com/tbrowne5/A-generalizable-Cas9-sgRNA-prediction-model-using-machine-transfer-learning)
#### [The original crisprHAL publication](https://doi.org/10.1038/s41467-023-41143-7)

The CRISPR/Cas9 nuclease from Streptococcus pyogenes (SpCas9) can be used with single guide RNAs (sgRNAs) as a sequence-specific antimicrobial agent and as a genome-engineering tool. However, current bacterial sgRNA activity models poorly predict SpCas9/sgRNA activity and are not generalizable, possibly because the underlying datasets used to train the models do not accurately measure SpCas9/sgRNA cleavage activity and cannot distinguish cleavage activity from toxicity. We solved this problem by using a two-plasmid positive selection system to generate high-quality biologically-relevant data that more accurately reports on SpCas9/sgRNA cleavage activity and that separates activity from toxicity. We developed a new machine learning architecture (crisprHAL) that can be trained on existing datasets and that shows marked improvements in sgRNA activity prediction accuracy when transfer learning is used with small amounts of high-quality data. The crisprHAL model recapitulates known SpCas9/sgRNA-target DNA interactions and provides a pathway to a generalizable sgRNA bacterial activity prediction tool.

## QUICK START

If you wish to run the model on your own nucleotide sequence follow parts 0 to 3.

If you wish to validate the listed models, follow parts 4 to 5.

## Sections of this guide:

Setting up and running the model to predict sgRNA activities:
* 0: Model requirements ```Time: 1-10 minutes```
* 1: Running the model test ```Runtime: 10 seconds```
* 2: Processing nucleotide sequences into model inputs ```Runtime: 1 second```
* 3: Predicting with the model ```Runtime: 10 seconds```

Additional information and methods: 
* 4: Preparing your own model input files & comparing predictions
* 5: Validating the trained models ```Variable runtime```
* 6: Data availability and processing
* 7: Citations
* 8: How to cite this model

## 0: Requirements

These are in a file called requirements.txt and should be in the working directory.
```
python>=3.8.8
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

These can be instantiated within a conda environment: ```Time: 1-10 minutes```

```
conda create --name HAL
conda activate HAL
conda install --file requirements.txt
```

This installation has been tested in Ubuntu 20.04.4 and Mac OSX 10.14.5, but has not been tested on Windows.

## 1: Run model test
```
python crisprHAL.py
```
Test our TevSpCas9 model with an example SpCas9 dataset of 7821 unique sgRNA target sites from Guo et al. 2018. 

Success here is that the model runs without error, showing that it is installed correctly. ```Runtime: ~10 seconds```



## 2: Process a fasta file of nucleotide sequence(s) into sgRNA target model inputs

This will take an input nucleotide fasta file and identifies potential sgRNA sequences for evaluation. The output will be
a .csv file containing the predicted guides. This can be used as input for the prediction step. ```Runtime: ~1 second```

* **Input Nucleotide File**: One single-line or multi-line fasta-formatted nucleotide sequence starting with a ">IDENTIFIER"
* **Output**: 28 nucleotide sequences in a CSV file appropriate as an input to the model

Composition of the 28 nucleotide inputs:
* 20 nucleotide target site, ie: CTCGATTGAGGGGCTGGGAA
* 3 nucleotide NGG PAM, ie: TGG
* 5 nucleotides downstream, ie: GTGAT
* Total: CTCGATTGAGGGGCTGGGAATGGGTGAT

Example input file and run shown below with a phiX174 genome:
```
>Sequence1
TCGAGCATGCATCTAGAGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCC
...etc

# python process_fasta.py [Input Nucleotide File] [Output CSV File]
python process_fasta.py phiX174.fna phiX174_sgRNAs.csv

#output: phiX174_output.csv
```

## 3: Predict with model

This will take the file of the predicted sgRNA sequences from above and assign a score. Higher scores are better!
The output is a .csv file named OUTPUT_[inputfile] and contains the sgRNA and the score. ```Runtime: ~10 seconds [+2-4 seconds/10,000 sites]```

* **Enzyme**: "TevSpCas9" or "SpCas9"
* **Input**: CSV file input name; created in section 2 or matching the required format (Format: Section 4)
* **Optional Compare**: "Compare" if your CSV file contains scores for comparison to model predictions (Format: Section 4)
* **Output**: Tab-deliminated (TSV) file containing the 28 nucleotide sequence and predicted Cas9 activity value

Example run with the phiX174 predicted sgRNA seqeunces
```
# python crisprHAL.py [Enzyme] [Input file csv] [Optional compare]
python crisprHAL.py TevSpCas9  phiX174_sgRNAs.csv
# output: OUTPUT_phiX174_sgRNAs.csv
```


Example command with prediction only, no "compare" option:
```
python crisprHAL.py TevSpCas9 test_dataset.csv
```

Example command with prediction and "compare" option for prediction comparison:
```
python crisprHAL.py SpCas9 test_dataset.csv compare
```


## 4: Preparing your own input CSV Files

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

## 5: Validate the training of the model

This will assess whether the training model is working. It will not change the model used for predictions.

Perform 5-fold cross validation with the TevSpCas9 dataset transfer learning from the eSpCas9 base model: ```Runtime: ~20 seconds```
```
python crisprHAL.py train TevSpCas9
```

Perform 5-fold cross validation with the SpCas9 dataset transfer learning from the eSpCas9 base model: ```Runtime: ~20 seconds```
```
python crisprHAL.py train SpCas9
```

Perform 80:20 train-test split with the eSpCas9 dataset (Guo et al. 2018) used as the base model\*: ```Runtime: 1-10 minutes```
```
python crisprHAL.py train eSpCas9
```
\*An 80% training & 20% test split was used for base model generation, and therefore has been included in place of the 5-fold cross validation tests used for the TevSpCas9 and SpCas9 enzyme transfer learning-based models.

## 6: Data availability and processing

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699

NCBI SRA Bioproject: PRJNA939699

Information about data processing can be found under crisprHAL/data/processing.


## 7: Citations

* Guo, J. et al. Improved sgRNA design in bacteria via genome-wide activity profiling. Nucleic acids research **46**, 7052–7069 (2018).
* Abadi, M. et al. TensorFlow: Large-scale machine learning on heterogeneous systems (2015). Software available from tensorflow.org.
* Fernandes, A. D. et al. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome **2**, 1–13 (2014).
* Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in python. Nat. Methods **17**, 261–272 (2020).

## 8: How to cite crisprHAL

Ham, D.T., Browne, T.S., Banglorewala, P.N. et al. A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets. **Nat Commun** *14*, 5514 (2023). https://doi.org/10.1038/s41467-023-41143-7
