
import re
import numpy as np
import pandas as p
import tensorflow as tf
import tensorflow.keras as k
from tensorflow.keras.layers import Dropout, MaxPooling1D, Dropout
from sklearn.model_selection import train_test_split, KFold
from Bio.SeqUtils import MeltingTemp as mt  
from scipy.stats import spearmanr
from statistics import stdev
import os
import sys
import encoder as e

os.environ["CUDA_VISIBLE_DEVICES"] = "1" 
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

def get_melt_temp(data):
    return np.array([mt.Tm_NN(guide) for guide in data.index])

def prepare_data_for_input(data, encoder=None, encoded_name=None, seqpadup=10, seqpaddown=5):
    if encoder:
        data[encoded_name] = [encoder.encode("A"*seqpadup + guide + "A"*seqpaddown) for guide in data.index]
    return data

def formatinputs(dataframe_input):
    selective=dataframe_input
    nrows, ncols = selective.shape
    temps = get_melt_temp(selective)
    long_temps = np.repeat(temps, ncols)
    selective = prepare_data_for_input(selective, e.BinaryEncoder, 'Binary', seqpadup=10, seqpaddown=5)
    binary_guides = np.repeat(selective.pop('Binary'),ncols)
    xprior = np.array([np.array(guide, dtype=np.float64) for guide in binary_guides])
    y_train = selective.values.reshape((nrows*ncols),)
    X_train = np.zeros((len(xprior),len(xprior[0]),5))
    for i in range(0, len(xprior)):
        for j in range(0, len(X_train[0])): 
            X_train[i][j] = np.append(xprior[i][j], long_temps[i]/100)
    y_train = np.array(y_train.astype(np.float))    
    return X_train, y_train


def get_ttsplit(data,ttsplitstate=1,tst_size=0.2,return_pre_format_values=False,prediction=False):
    dataset = p.read_csv(data, header=None, skiprows=0,).dropna(how="all")
    if prediction: dataset['scores'] = 0
    nparray = dataset.to_numpy()
    nparray = nparray[1:]
    if tst_size==0:
        df_all = p.DataFrame(data=nparray[:,1:],index=nparray[:,0])
        X_all,y_all=formatinputs(df_all)
        if return_pre_format_values==False: return X_all,y_all
        else: return df_all,X_all,y_all
    else:
        All_train, All_test = train_test_split( nparray, test_size=tst_size, random_state=ttsplitstate)
        df_train = p.DataFrame(data=All_train[:,1:],index=All_train[:,0])
        df_test = p.DataFrame(data=All_test[:,1:],index=All_test[:,0])
        X_train,y_train=formatinputs(df_train)
        X_test,y_test=formatinputs(df_test)
        if return_pre_format_values==False: return X_train, X_test, y_train, y_test
        else: return df_train, df_test, X_train, X_test, y_train, y_test


def write_prediction(y_pred, df_test, predictionfilename):
    predictionwrite = open(predictionfilename,"w+")
    np_train = df_test.index.to_numpy()
    for i in range(0,len(np_train)):
        predictionwrite.write(np_train[i] + "\t" + str(y_pred[i]).replace(' [','').replace('[', '').replace(']', '') + "\n")
    predictionwrite.close()


def main(fileoutput="NULL", train=False, compare=False, model="Tev", inputdata="X", seqstart=10, seqend=38, drop_rate=0.3, CNN_filters=128, window_size=3, CNN_drop=0.3, conv1D_padding="same", CNN_dense1=128, CNN_dense2=64, maxpool1D_padding="same", RNN_size=128, RNN_dense1=128, RNN_dense2=64, CNN_RNN_drop=0.3):

    if train:
        if inputdata=="eSpCas9":
            i = k.Input(shape=(43,5), name="Input")

            # Seqstart and seqend variables used to identify optimal sequence length for model performance
            x = i[:,seqstart:seqend,0:4]

            # Multi-layer CNN
            c1 = k.layers.Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c1")(x)
            l1 = k.layers.LeakyReLU(name="l1")(c1)
            p1 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p1")(l1)
            dr1 = Dropout(CNN_drop, name="dr1")(p1)

            c2 = k.layers.Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c2")(dr1)
            l2 = k.layers.LeakyReLU(name="l2")(c2)
            p2 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p2")(l2)
            dr2 = Dropout(CNN_drop, name="dr2")(p2)

            c3 = k.layers.Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c3")(dr2)
            l3 = k.layers.LeakyReLU(name="l3")(c3)
            p3 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p3")(l3)
            dr3 = Dropout(CNN_drop, name="dr3")(p3)
            
            c4 = k.layers.Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c4")(dr3)
            l4 = k.layers.LeakyReLU(name="l4")(c4)
            p4 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p4")(l4)
            dr4 = Dropout(CNN_drop, name="dr4")(p4)

            f = k.layers.Flatten(name="f")(dr4)

            d1 = k.layers.Dense(CNN_dense1, name="d1")(f)
            ld1 = k.layers.LeakyReLU(name="ld1")(d1)
            drd1 = Dropout(CNN_drop, name="drd1")(ld1)
            d2 = k.layers.Dense(CNN_dense2, name="d2")(drd1)
            ld2 = k.layers.LeakyReLU(name="ld2")(d2)
            drd2 = Dropout(CNN_drop, name="drd2")(ld2)
            x_1d_o = k.layers.Dense(1, name="x_1d_o")(drd2)

            # Bidirectional Gated Recurrent Unit Branch
            r1 = k.layers.Bidirectional(k.layers.GRU(RNN_size, kernel_initializer='he_normal', dropout=drop_rate, recurrent_dropout=0.2), name="r1")(dr1)
            
            f_LSTM = k.layers.Flatten(name="f_LSTM")(r1)
            
            d1_LSTM = k.layers.Dense(RNN_dense1, name="d1_LSTM")(f_LSTM)
            ld1_LSTM = k.layers.LeakyReLU(name="ld1_LSTM")(d1_LSTM)
            drd1_LSTM = Dropout(CNN_RNN_drop, name="drd1_LSTM")(ld1_LSTM)
            d2_LSTM = k.layers.Dense(RNN_dense2, name="d2_LSTM")(drd1_LSTM)
            ld2_LSTM = k.layers.LeakyReLU(name="ld2_LSTM")(d2_LSTM)
            drd2_LSTM = Dropout(CNN_RNN_drop, name="drd2_LSTM")(ld2_LSTM)
            x_LSTM_o = k.layers.Dense(1, name="x_LSTM_o")(drd2_LSTM)
            
            # Concatenated Output Layers
            o_c = k.layers.concatenate([x_1d_o,x_LSTM_o], axis=1, name="o_c")
            o = k.layers.Dense(1, activation="linear", name="o")(o_c)

            # Model Initiation
            m = k.Model(inputs=i, outputs=o)
            compile_options={"optimizer": "adam", "loss": "mean_squared_error"}
            m.compile(**compile_options)

            # Training the eSpCas9 data model from which the TevSpCas9 & SpCas9 models transfer learn
            dataX_train, dataX_test, datay_train, datay_test = get_ttsplit("data/eSpCas9.csv",1,0.2)
            m.fit(dataX_train,datay_train,epochs=40,batch_size=200,verbose=1)
            result = spearmanr(m.predict(dataX_test),datay_test,axis=0)
            print("\nSpearman ranked correlation coefficient: " + str(result[0]))
            
        else:
            # Test the model with TevSpCas9 data under 5-fold cross validation
            if "Tev" in indata or "tev" in indata:
                dataX, datay = get_ttsplit("data/TevSpCas9.csv",1,0)
            # Test the model with SpCas9 data under 5-fold cross validation
            elif "SpCas9" in indata or "Cas9" in indata or "cas9" in indata:
                dataX, datay = get_ttsplit("data/SpCas9.csv",1,0)
            else:
                print("ERROR: Invalid training option, please choose one of: TevSpCas9, SpCas9, or eSpCas9.")
            print("Testing the " + indata + " model.")

            kf = KFold(n_splits=5, shuffle=True, random_state=1)
            results=[]
            compile_options={"optimizer": "adam", "loss": "mean_squared_error"}

            # Running the 5-fold cross validation procedure
            for train_index, test_index in kf.split(dataX):
                cvX_train, cvX_test = dataX[train_index], dataX[test_index]
                cvy_train, cvy_test = datay[train_index], datay[test_index]
                m = k.models.load_model("h5/eSpCas9.h5")
                for i in range(0,len(m.layers)):
                    if i not in [18, 22, 28, 34, 35]:
                        m.layers[i].trainable = False
                m.compile(**compile_options)
                m.fit(cvX_train, cvy_train, epochs=4, batch_size=20, verbose=0)
                result = spearmanr(m.predict(cvX_test, verbose=0),cvy_test)
                results.append(result[0])
                #print(result)
            print("\nMean 5-fold cross validation score: " + str(sum(results)/5))
            print("Standard deviation of scores: " + str(stdev(results)))
            print("\nMeasurement performed by Spearman ranked correlation")

    else:
        # Load either the TevSpCas9 or SpCas9 model
        if "Tev" in model or "tev" in model:
            m = k.models.load_model("h5/TevSpCas9.h5")
        elif ("SpCas9" in model or "Cas9" in model or "cas9" in model) and "Tev" not in model and "eSpCas9" not in model:
            m = k.models.load_model("h5/SpCas9.h5")
        else:
            print("ERROR: No correct model specified, please choose either: TevSpCas9 or SpCas9")
            return
        
        # Option: Comparison of the model predictions to prior scores
        if compare:
            input_data, input_Xall, input_yall = get_ttsplit(inputdata,1,0,True)
            input_pred = m.predict(input_Xall, verbose=0)
            spearman_corr = spearmanr(input_pred,input_yall,axis=0)
            print("Spearman ranked correlation coefficient: " + str(spearman_corr[0]))
        else:
            input_data, input_Xall, input_yall = get_ttsplit(inputdata,1,0,True,prediction=True)
            input_pred = m.predict(input_Xall, verbose=0)
        if fileoutput != "NULL":
            write_prediction(input_pred, input_data, fileoutput)

compare=False
train=False
if len(sys.argv) > 1:
    modelname=str(sys.argv[1])
    if modelname=="train" or modelname=="Train": train=True
    if len(sys.argv) > 2: indata=str(sys.argv[2])
    if len(sys.argv) > 3 and ("compare" in str(sys.argv[3]) or "Compare" in str(sys.argv[3])): compare=True
    outfile="Output_" + indata
    main(inputdata=indata,model=modelname,fileoutput=outfile,compare=compare,train=train)
else:
    print("Beginning the crisprHAL.py model test")
    print("\nRunning the TevSpCas9 model with an example SpCas9 dataset of 7821 sgRNAs from Guo et al. 2018")
    main(inputdata="test_dataset.csv",model="TevSpCas9",fileoutput="Output_TevSpCas9_test_dataset.csv",compare=True,train=False)
    print("\nRunning the SpCas9 model with an example SpCas9 dataset of 7821 sgRNAs from Guo et al. 2018")
    main(inputdata="test_dataset.csv",model="SpCas9",fileoutput="Output_SpCas9_test_dataset.csv",compare=True,train=False)
