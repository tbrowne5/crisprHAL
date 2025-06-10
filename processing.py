
import numpy as np
from scipy.stats import spearmanr, pearsonr
from models import modelVersionInputLength, modelVersionTrainingData, modelVersionTestingData

class processing:
    
    def __init__(self):
        pass

    def onehotencode(self, nucleotideEncoding):
        
        nt_to_ohe_dict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "M": [0, 0, 0, 0], "N": [0, 0, 0, 0]}
        onehotencoding = []
        for sequence in nucleotideEncoding: onehotencoding.append(np.array([nt_to_ohe_dict[nt.upper()] for nt in sequence]))
        return np.array(onehotencoding)

    def find_targets(self, sequence, inputParameters, circular=True, reverseComplement=True):
        
        if len(sequence) < inputParameters[0]:
            print(f"Error: Input sequence '{sequence}' is shorter than the required length of {inputParameters[0]} bases.")
            return []
        else:
            sequence = sequence.upper().replace("U", "T").replace(" ", "").replace("\n", "")
            if circular: sequence = sequence + sequence[:inputParameters[0] - 1]
            targets = []
            # Find NGG targets
            for i in range(20+inputParameters[1], len(sequence) - inputParameters[2] + 1):
                if sequence[i+1:i+3] == "GG":
                    target = sequence[i-(20+inputParameters[1]):i+inputParameters[2]]
                    if len(target) == inputParameters[0]:
                        targets.append(target)
                    else:
                        print(len(target))
                        break
            # Find NGG targets in reverse complement
            if reverseComplement:
                sequence = sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))
                for i in range(20+inputParameters[1], len(sequence) - inputParameters[2] + 1):
                    if sequence[i+1:i+3] == "GG":
                        target = sequence[i-(20+inputParameters[1]):i+inputParameters[2]]
                        if len(target) == inputParameters[0]:
                            targets.append(target)
                        else:
                            print(len(target))
                            break
            return targets

    def process_fasta(self, fastaFile, inputParameters=[37, 3, 14], circular=True, reverseComplement=True):
        
        # Read the FASTA file
        inputSequences = []
        with open(fastaFile, 'r') as file:
            sequence = ""
            for line in file:
                if line.startswith('>'):
                    if len(sequence) > 0:
                        inputSequences += self.find_targets(sequence, inputParameters, circular, reverseComplement)
                        # ALSO NEED TO ADD IN THE REVERSE COMPLEMENT OF THE SEQUENCE
                        sequence = ""
                else:
                    sequence += line.strip()
            inputSequences += self.find_targets(sequence, inputParameters, circular, reverseComplement)
        
        print(f"Found {len(inputSequences)} target sequences in the FASTA file.")
        return inputSequences, self.onehotencode(inputSequences), None

    def process_csv(self, csvFile, predictionColumn=False):
        
        # Check if comma or tab separated via file name
        if csvFile.endswith('.csv'): delimiter = ','
        elif csvFile.endswith('.tsv'): delimiter = '\t'
        else: raise ValueError("Unsupported file format. Please provide a .csv or .tsv file.")
        
        # Read the CSV file
        inputSequences = []
        scores = []
        
        with open(csvFile, 'r') as file:
            for line in file:
                if line.strip():
                    targetAndScore = line.strip().split(delimiter)
                    inputSequences.append(targetAndScore[0])
                    if predictionColumn and len(targetAndScore) > 1: scores.append(float(targetAndScore[1]))

        oneHotEncoded = self.onehotencode(inputSequences)

        if predictionColumn: return inputSequences, oneHotEncoded, np.array(scores)
        else: return inputSequences, oneHotEncoded, None

    def read_training_data(self, modelName):
        return self.process_csv(modelVersionTrainingData[modelName], predictionColumn=True)

    def read_testing_data(self, modelName):
        return self.process_csv(modelVersionTestingData[modelName], predictionColumn=True)

    def read_input(self, modelName, inputFile, compare, circular=True):

        # If no parameters specified, the program will provide hold-out testing data performance
        if inputFile == None:
            return self.process_csv(modelVersionTestingData[modelName], predictionColumn=True)
        # Recommended input format for model use (no fiddling with input sequence lengths)
        elif inputFile.endswith('.fasta') or inputFile.endswith('.fa'):
            return self.process_fasta(inputFile, modelVersionInputLength[modelName], circular)
        # Method for validating model performance, requires explicit input sequence lengths
        elif inputFile.endswith('.csv') or inputFile.endswith('.tsv'):
            return self.process_csv(inputFile, predictionColumn=compare)
        else:
            raise ValueError("Unsupported file format, please provide an input of the format: .fasta, .fa, .csv, or .tsv")
    
    def compare_predictions(self, predictions, actual, returnScores=False, message=None):
        
        if len(predictions) != len(actual):
            raise ValueError("Predictions and actual scores must have the same length.")
        
        # Calculate Spearman and Pearson correlations
        predictions = np.array(predictions).flatten()
        spearman_corr, _ = spearmanr(predictions, actual, axis=0)
        pearson_corr, _ = pearsonr(predictions, actual, axis=0)

        if message is not None: print(message)

        print(f"Spearman correlation: {spearman_corr:.4f}")
        print(f"Pearson correlation: {pearson_corr:.4f}")

        if returnScores: return spearman_corr, pearson_corr
    
    def write_predictions(self, inputSequences, predictions, outputFile=None, inputFile=None, inputScores=None):

        # If output file is not specified, use inputfile name with "_predictions.csv" suffix
        # If neither inputFile nor outputFile is specified, write to "example_predictions.csv"
        if outputFile is None:
            if inputFile is not None:
                outputFile = inputFile.rsplit('.', 1)[0] + "_predictions.csv"
            else:
                outputFile = "example_predictions.csv"
        # Predictions and input scores are numpy arrays
        predictions = np.array(predictions).reshape(-1, 1)
        if inputScores is not None:
            inputScores = [str(score[0]) for score in np.array(inputScores).reshape(-1, 1)]
        # Write predictions to a CSV file
        with open(outputFile, 'w') as file:
            if inputScores is not None:
                file.write("sgRNA,Score,Predictions\n")
                for seq, score, pred in zip(inputSequences, inputScores, predictions):
                    file.write(f"{seq},{score},{pred[0]}\n")
            else:
                file.write("sgRNA,Predictions\n")
                for seq, pred in zip(inputSequences, predictions):
                    file.write(f"{seq},{pred[0]}\n")

