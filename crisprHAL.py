
import re
import numpy as n
import numpy as np
import pandas as p
import tensorflow as tf
import tensorflow.keras as k
from tensorflow.keras.layers import InputLayer
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import LeakyReLU
from tensorflow.keras.layers import Dropout
from tensorflow.keras.models import load_model
from tensorflow.keras.layers import MaxPooling1D, Dropout, Activation
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier, KerasRegressor
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, GridSearchCV, KFold
from sklearn import datasets, svm
from Bio.SeqUtils import MeltingTemp as mt  
from scipy.stats import gmean  # I added this
from scipy.stats import spearmanr
from statistics import stdev
import array as arr
import os

import sys

import encoder as e

os.environ["CUDA_VISIBLE_DEVICES"] = "1" 
transfer_learning_model=0

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


def get_ttsplit (data,ttsplitstate=1,tst_size=0.2,return_pre_format_values=False,prediction=False):
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


def main(fileoutput="NULL", compare=False, model="Tev", seqstart=10, seqend=38, basedata='X', transferdata='X', inputdata='X', epochs=40, batch_size=200, tl_epochs=4, drop_rate=0.3, CNN_filters=128, window_size=3, CNN_drop=0.3, conv1D_padding="same", CNN_dense1=128, CNN_dense2=64, maxpool1D_padding="same", RNN_size=128, RNN_dense1=128, RNN_dense2=64, CNN_RNN_drop=0.3, verbose=0):
    i = k.Input(shape=(43,5), name="Input")

    # Using seqstart and seqend inputs to identify optimal sequence length for model performance
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
    r1 = k.layers.Bidirectional(k.layers.GRU(RNN_size, kernel_initializer='he_normal', dropout=drop_rate, recurrent_dropout=0.2), name="r1")(dr1) #dr1_LSTM) #LSTM(64,activation=None,dropout=0.3,recurrent_dropout=0.2)(dr1_LSTM) #dr1_LSTM)
    
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

    if "Tev" in model or "tev" in model:
        m = k.models.load_model("/Volumes/bin/CrisprHAL_prediction_h5s/CrisprHAL_predictions/tlrn_mergestart_Tev_std_72_1_1_CrisprHAL_stdev_eSpTev_HH_128_128_128_64_branch_10_to_38_eSp_stdev_scaled_43_UniqCas9_43_unscaled_difference_Tev_stdev_scaled_43_4_allData.h5")
    elif ("SpCas9" in model or "Cas9" in model or "cas9" in model) and "Tev" not in model and "eSpCas9" not in model:
        m = k.models.load_model("/Volumes/bin/CrisprHAL_prediction_h5s/CrisprHAL_predictions/tlrn_mergestart_pTCas_s_56_1_1_CrisprHAL_stdev_eSppTCas_F_128_128_128_64_branch_10_to_38_eSp_stdev_scaled_43_UniqCas9_43_unscaled_difference_pTCas_stdev_scaled_43_6_allData.h5")
    else:
        print("No correct model specified, please choose either: TevSpCas9 or SpCas9")
        return
    
    if compare:
        input_data, input_Xall, input_yall = get_ttsplit(inputdata,1,0,True)
        input_pred = m.predict(input_Xall, verbose=0)
        spearman_corr = spearmanr(input_pred,input_yall,axis=0)
        print(spearman_corr[0])
    else:
        input_data, input_Xall, input_yall = get_ttsplit(inputdata,1,0,True,prediction=True)
        input_pred = m.predict(input_Xall, verbose=0)
    if fileoutput != "NULL":
        write_prediction(input_pred, input_data, fileoutput)

compare=False
if len(sys.argv) > 1: modelname=str(sys.argv[1])
else: modelname="Tev"

if len(sys.argv) > 2: indata=str(sys.argv[2])
else:
    print("Running model with an example SpCas9 dataset of 7821 sgRNAs from Guo et al. 2018")
    indata="test_dataset.csv"
    compare=True

if len(sys.argv) > 3:
    if "compare" in str(sys.argv[3]) or "Compare" in str(sys.argv[3]): compare=True
outfile="Output_" + indata

main(inputdata=indata,model=modelname,fileoutput=outfile,compare=compare)
