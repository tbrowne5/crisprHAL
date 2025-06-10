# crisprHAL: A suite of bacterial sgRNA/Cas9 activity prediction models
The up-to-date models for bacterial SpCas9, TevSpCas9, eSpCas9, and SaCas9 nuclease activity prediction.

## Links to the crisprHAL repositories, papers, and web tool:
* [Up-to-date crisprHAL prediction tool for use](https://github.com/tbrowne5/crisprHAL) — **YOU ARE HERE**
* [Online crisprHAL prediction tool](https://crisprhal.streamlit.app/)
* [crisprHAL 2.0 paper repository](https://github.com/tbrowne5/Better-data-for-better-predictions-data-curation-improves-deep-learning-for-sgRNA-Cas9-prediction/)
* crisprHAL 2.0 pre-print (Available soon)
* [crisprHAL SaCas9 paper repository](https://github.com/tbrowne5/Adenine-methylated-PAM-sequences-inhibit-SaCas9-activity)
* crisprHAL SaCas9 pre-print (Available soon)
* [crisprHAL SpCas9 paper repository](https://github.com/tbrowne5/A-generalizable-Cas9-sgRNA-prediction-model-using-machine-transfer-learning)
* [crisprHAL SpCas9 publication](https://doi.org/10.1038/s41467-023-41143-7)

## Available Prediction Models:

* SpCas9: ```-m SpCas9``` for SpCas9 and TevSpCas9 activity prediction
* eSpCas9: ```-m eSpCas9``` for eSpCas9 activity prediction
* SaCas9: ```-m TevSpCas9``` for SaCas9 and TevSaCas9 activity prediction
* WT-SpCas9: ```-m WT-SpCas9``` a secondary SpCas9 prediction model

## QUICK START

If you wish to run the model on your own nucleotide sequence follow parts 0 to 3.

If you wish to validate the training of a model, follow parts 4 to 5.

If you wish to simply obtain predictions, you can do so easily through the [crisprHAL website](https://crisprHAL.streamlit.app).

As this repository contains the most up-to-date models, if you wish to test a model from a specific paper, please use the paper-specific repository link listed above.

## SECTIONS OF THIS GUIDE:

Setting up and running the model to predict sgRNA/Cas9 activities:
* 0: Model requirements ```Time: 1-10 minutes```
* 1: Running the model test ```Runtime: 10 seconds```
* 2: Understanding available model options
* 3: Predicting with the model ```Runtime: 1-10 seconds```

Additional information and methods: 
* 4: Preparing your own model input files & comparing predictions
* 5: Validating the trained models ```Variable runtime```
* 6: Data availability and processing
* 7: Citations
* 8: How to cite this model

## 0: Requirements

These are in a file called requirements.txt and should be in the working directory.
```
python
tensorflow
```

These can be instantiated within a conda environment: ```Time: 1-10 minutes```

```
conda create --name HAL
conda activate HAL
conda install --file requirements.txt
```

This installation has been tested in Ubuntu 20.04.4 and Mac OSX 10.14.5, but has not been tested on Windows.


## 1: Run the model test
```
python crisprHAL.py
```
Test our primary SpCas9/TevSpCas9 prediction model: crisprHAL Tev

Success here is that the model runs without error, showing that it is installed correctly. ```Runtime: ~10 seconds```


## 2: Understand options

```
python crisprHAL.py [options]
```
```
--model, -m   [TevSpCas9, eSpCas9, WT-SpCas9]
```
Model (default=TevSpCas9): specify the model to be used. TevSpCas9 should be used for both TevSpCas9 and SpCas9 predictions. WT-SpCas9 should only be used for crisprHAL WT validation.
```
--input, -i   [Input file path]
```
Input: crisprHAL accepts three types of input files: FASTA (.fasta and .fa), CSV (.csv), and TSV (.tsv). If no input is specified, the model will default to testing on its respective hold-out set.

```
--output, -o  [Output file path]
```
Output: specify the output path and file name of choice. If no output is specified, output file will have the input file name with the prefix: "_predictions.txt"
```
--circular
```
Circular (default=FALSE): specific to FASTA inputs; specifies that the input sequence should be treated as circular DNA rather than linear.
```
--compare, -c
```
Compare (default=FALSE): specific to CSV/TSV inputs; specifies that the input file contains a second column with scores for comparison. Outputs Spearman and Pearson correlation coefficients between predictions and provided scores, and writes both the predictions and scores to the output file.
```
--train, -t
```
Train (default=FALSE): train the model specified by ```--model/-m [modelName]``` (default=TevSpCas9) using the corresponding training dataset in the data/ directory.
```
--epochs, -e  [Integer epoch value]
```
Epochs: specify the number of training epochs to be run when training the model. By default each model uses its respective 5CV-identified epochs. Without the ```--train/-t``` flag, ```--epochs/-e``` will do nothing.
```
--help, -h
```
Help: prints available options.


## 3: Predict with the model

To run a test of the model's predictions, please run the command:
```
python crisprHAL.py -m TevSpCas9 -i sample.fa
```
This will: 1) Identify all sgRNA targets in the sample.fa file, 2) Predict Cas9 activity for each target, and 3) Write the targets and predicted activities to a file called sample_predictions.txt

The full list of options for FASTA input-based model predictions are:
```
python crisprHAL.py --model [TevSpCas9/eSpCas9/WT-SpCas9] --input [input file path] --output [optional output file path] --circular
```


## 4: Preparing your own input CSV Files

Input CSV file for prediction only, no "Compare" option:
```
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
