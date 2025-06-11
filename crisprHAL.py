import sys
import os
from models import models, modelVersionDefaultEpochs
from processing import processing

training = False
modelName = "TEVSPCAS9"
modelNames = ["TEVSPCAS9","TEVCAS9", "TEV", "SPCAS9", "WTSPCAS9", "WT-SPCAS9", "WILD-TYPE-SPCAS9","WILD-TYPE","WILDTYPE", "WT", "ESPCAS9", "ESP","TEVSACAS9","SACAS9","SACAS","SA"]
epochs = None
inputFile = None
outputFile = None
compare = False
circularInput = False
summary = False

# PERHAPS MOVE THIS TO PROCESSING AND HAVE A VERY SIMPLE CRISPRHAL.PY FILE!!!

def parse_args(args):
    global training, modelName, modelNames, epochs, inputFile, outputFile, compare, circularInput, summary

    for i in range(len(args)):
        if args[i] == "--train" or args[i] == "-t": training = True
        elif args[i] == "--input" or args[i] == "-i": inputFile = args[i + 1]
        elif args[i] == "--output" or args[i] == "-o": outputFile = args[i + 1]
        elif args[i] == "--circular": circularInput = True
        elif args[i] == "--compare" or args[i] == "-c": compare = True
        elif args[i] == "--epochs" or args[i] == "-e": epochs = int(args[i + 1])
        elif args[i] == "--summary" or args[i] == "-s": summary = True
        elif args[i] == "--model" or args[i] == "-m" or args[i] == "--enzyme":
            if args[i + 1].upper() in modelNames:
                if args[i + 1].upper() in ["TEVSPCAS9","TEVCAS9", "TEV", "SPCAS9"]: modelName = "TEVSPCAS9"
                if args[i + 1].upper() in ["WT-SPCAS9", "WILD-TYPE-SPCAS9","WILD-TYPE","WILDTYPE", "WT"]: modelName = "WT-SPCAS9"
                if args[i + 1].upper() in ["ESPCAS9", "ESP"]: modelName = "ESPCAS9"
                if args[i + 1].upper() in ["TEVSACAS9","SACAS9","SACAS","SA"]: modelName = "TEVSACAS9"
            else:
                print(f"Error: Model '{args[i + 1]}' is not recognized. Available models: TevSpCas9, eSpCas9, and WT-SpCas9.\nFor up-to-date and other models including SaCas9 please use: github.com/tbrowne5/crisprHAL")
                sys.exit(1)
        elif args[i] == "--help" or args[i] == "-h":
            print("Thank you for using crisprHAL. This program offers bacterial sgRNA/Cas9 activity predictions for several nucleases. For typical use please provide a fasta file input (.fasta/.fa).")
            print("\nThe nuclease options are:\n  • SpCas9: Streptococcus pyogenes Cas9 with NGG PAM\n  • TevSpCas9: I-TevI linker SpCas9 dual nuclease with NGG PAM\n  • eSpCas9: enhanced Specificity SpCas9 with NGG PAM\n  • SaCas9: Staphylococcus aureus Cas9 with NNGRRN PAM\n  • TevSaCas9: I-TevI linker SaCas9 dual nuclease with NGG PAM\n")
            print("Usage: python crisprHAL.py [options]")
            print("Options:")
            print("  --model, -m   [Nuclease]                Specify the model name (default: TevSpCas9)")
            print("  --input, -i   [Input file path]         Input file for prediction (fasta, csv, or tsv)")
            print("  --output, -o  [Output file path]        Output file for prediction results")
            print("  --circular                              Process fasta as a circular input sequence")
            print("  --compare, -c                           Compare predictions with scores in the input file second column")
            print("  --train, -t                             Train the model specified")
            print("  --epochs, -e  [Integer epoch value]     Specify number of epochs for training (default: model-specific)")
            print("  --summary, -s                           Print the model architecture summary")
            print("  --help, -h                              Show this help message")
            sys.exit(0)
    
    if training == True and inputFile is not None:
        print("Error: Please specify either training mode or an input file for prediction generation.")
        sys.exit(1)

# Read in arguments
# If no input specified, use default hold-out test set
# If no output specified, but input specified, strip file type annd add _predictions.csv

def run_model():
    global training, modelName, inputFile, outputFile, compare, epochs, circularInput, summary

    process = processing()
    model = models(modelName, summary)
    
    if training:
        print("Training model")
        trainingData = process.read_training_data(modelName)
        testingData = process.read_testing_data(modelName)
        if epochs is None: epochs = modelVersionDefaultEpochs[modelName]
        for i in range(0,epochs):
            model.train(trainingData[1], trainingData[2], epochs=1, batch_size=1024, verbose=1)
            process.compare_predictions(model.predict(testingData[1]), testingData[2], message=f"Epoch {i+1} on hold-out test set for {modelName}:")
    else:
        # Model name provides input sequence length for processing
        # If inputFile default of "None" is passed, the hold-out test set will be used instead
        # The compare flag indicates that the input file contains a second column of scores to be used for comparison
        if inputFile is None: compare = True
        inputSequences, encodedInputSequences, inputScores = process.read_input(modelName, inputFile, compare, circularInput)
        print(len(inputSequences))
        model.load_model(modelName)
        predictions= model.predict(encodedInputSequences)
        if compare:
            # Compare predictions with the second column of scores in the input file
            process.compare_predictions(predictions, inputScores)
        process.write_predictions(inputSequences, predictions, outputFile, inputFile, inputScores)

if __name__ == "__main__":
    parse_args(sys.argv[1:])
    run_model()
