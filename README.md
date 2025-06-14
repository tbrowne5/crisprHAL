# crisprHAL: A suite of bacterial sgRNA/Cas9 activity prediction models
The up-to-date models for bacterial SpCas9, TevSpCas9, eSpCas9, and SaCas9 nuclease activity prediction.

## Links to the crisprHAL repositories, papers, and web tool:
* [Up-to-date crisprHAL prediction tool for use](https://github.com/tbrowne5/crisprHAL) — **YOU ARE HERE**
* [Online crisprHAL prediction tool](https://crisprhal.streamlit.app/) ([Repository](https://github.com/tbrowne5/crisprHAL_streamlit))
* [crisprHAL 2.0 paper repository](https://github.com/tbrowne5/Better-data-for-better-predictions-data-curation-improves-deep-learning-for-sgRNA-Cas9-prediction/)
* crisprHAL 2.0 pre-print (Available soon)
* [crisprHAL SaCas9 paper repository](https://github.com/tbrowne5/Adenine-methylated-PAM-sequences-inhibit-SaCas9-activity)
* crisprHAL SaCas9 pre-print (Available soon)
* [crisprHAL SpCas9 paper repository](https://github.com/tbrowne5/A-generalizable-Cas9-sgRNA-prediction-model-using-machine-transfer-learning)
* [crisprHAL SpCas9 publication](https://doi.org/10.1038/s41467-023-41143-7)

## Available Prediction Models:

* SpCas9: ```-m SpCas9``` for SpCas9 and TevSpCas9 activity prediction
* eSpCas9: ```-m eSpCas9``` for eSpCas9 activity prediction
* SaCas9: ```-m TevSpCas9``` for SaCas9 and TevSaCas9 activity prediction (**In-Progress**)
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

Python and Tensorflow versions used:
```
python==3.12
tensorflow==2.19
```
Create a virtual environment to install Tensorflow and SciPy:
```
python3.12 -m venv HAL
source HAL/bin/activate
pip install tensorflow[and-CUDA]==2.19
pip install scipy==1.15
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
Options:
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

Preparing your own input CSV (or TSV) file requires more work, as the input must be the exact length required by the specific model. Different models have different target site adjacent sequence inclusion requirements as follows:

* SpCas9/TevSpCas9 prediction: ```37 nucleotides = 3 upstream, 20 sgRNA, 14 downstream (NGG PAM plus 11)```
* eSpCas9 prediction: ```406 nucleotides = 193 upstream, 20 sgRNA, 193 downstream (NGG PAM plus 190)```
* WT-SpCas9 prediction: ```378 nucleotides = 189 upstream, 20 sgRNA, 169 downstream (NGG PAM plus 166)```
* SaCas9/TevSaCas9 prediction: ```29 nucleotides = 21 sgRNA, 8 downstream (NNGRRN PAM plus 2)```

Input CSV file for prediction only, SpCas9/TevSpCas9 model input, no "Compare" option:
```
GCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTA
GATACCAGGATCTTGCCATCCTATGGAACTGCCTCGG
GTAACCATGCATCATCAGGAGTACGGATAAAATGCTT
CTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGG
AATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGG
GTGACGACTGAATCCGGTGAGAATGGCAAAAGCTTAT
AAAACGCCATTAACCTGATGTTCTGGGGAATATAAGG
CTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTT
```

Predict using the command: ```python crisprHAL.py -i [fileName.csv]``` with optional specification of model ```-m [modelName]``` and output path ```-o [outputPath]```

Input CSV file for prediction and prediction vs score comparison in column 2, SpCas9/TevSpCas9 model input:
```
GCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTA,1.17419584308754
GATACCAGGATCTTGCCATCCTATGGAACTGCCTCGG,-0.482944213762817
GTAACCATGCATCATCAGGAGTACGGATAAAATGCTT,0.421119225849831
CTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGG,0.654659434852531
AATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGG,-0.995596435390042
GTGACGACTGAATCCGGTGAGAATGGCAAAAGCTTAT,-2.07685671663022
AAAACGCCATTAACCTGATGTTCTGGGGAATATAAGG,-0.294359200889779
CTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTT,0.92532668013273
```

Predict using the command: ```python crisprHAL.py -i [fileName.csv] -c``` with optional specification of model ```-m [modelName]``` and output path ```-o [outputPath]```


## 5: Validate the training of the model

This will assess whether the training model is working. It will not change the model used for predictions.

Train the model for the model-specific default number of epochs. If no model is specified, crisprHAL Tev (SpCas9/TevSpCas9) will be trained.
```
python crisprHAL.py -t -m [modelName]
```

Traing the model for a specific number of epochs:
```
python crisprHAL.py -t -m [modelName] -e [integerValueForTrainingEpochs]
```

Example: train crisprHAL eSp (eSpCas9) for 50 epochs:
```
python crisprHAL.py -m eSpCas9 -t -e 50
```


## 6: Data availability and processing

Raw sequence reads from which the datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699 https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1260991 and https://www.ncbi.nlm.nih.gov/bioproject/PRJNA450978

NCBI SRA Bioproject: PRJNA939699

NCBI SRA BioProject: PRJNA1260991

NCBI SRA BioProject: PRJNA450978


## 7: Citations

* Guo, J. et al. Improved sgRNA design in bacteria via genome-wide activity profiling. Nucleic acids research **46**, 7052–7069 (2018).
* Abadi, M. et al. TensorFlow: Large-scale machine learning on heterogeneous systems (2015). Software available from tensorflow.org.
* Fernandes, A. D. et al. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome **2**, 1–13 (2014).
* Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in python. Nat. Methods **17**, 261–272 (2020).

## 8: How to cite crisprHAL

Ham, D.T., Browne, T.S., Banglorewala, P.N. et al. A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets. **Nat Commun** *14*, 5514 (2023). https://doi.org/10.1038/s41467-023-41143-7
