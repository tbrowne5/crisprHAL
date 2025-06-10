
import tensorflow.keras as k
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Input, Dense, Conv1D, MaxPooling1D, Conv2D, MaxPooling2D, Dropout, Bidirectional, GRU, LeakyReLU, Flatten, concatenate, BatchNormalization
from tensorflow.keras.models import Model
#from keras_gradient_noise import add_gradient_noise

modelVersionInputLength = {
    "TEVSPCAS9": [37, 3, 14],
    "ESPCAS9": [406, 193, 193],
    "WT-SPCAS9": [378, 189, 169],
    "TEVSACAS9": [29, 1, 8]
    }

modelVersionPath = {
    "TEVSPCAS9": "models/TevSpCas9.keras",
    "ESPCAS9": "models/eSpCas9.keras",
    "WT-SPCAS9": "models/WT-SpCas9.keras",
    "TEVSACAS9": "models/TEVSACAS9.h5"
    }

modelVersionTrainingData = {
    "TEVSPCAS9": "data/TevSpCas9_training_data.csv",
    "ESPCAS9": "data/eSpCas9_training_data.csv",
    "WT-SPCAS9": "data/WT-SpCas9_training_data.csv",
    "TEVSACAS9": "data/TevSaCas9_training_data.csv"
    }

modelVersionTestingData = {
    "TEVSPCAS9": "data/TevSpCas9_testing_data.csv",
    "ESPCAS9": "data/eSpCas9_testing_data.csv",
    "WT-SPCAS9": "data/WT-SpCas9_testing_data.csv",
    "TEVSACAS9": "data/TevSaCas9_testing_data.csv"
    }

modelVersionDefaultEpochs = {
    "TEVSPCAS9": 85,
    "ESPCAS9": 74,
    "WT-SPCAS9": 48,
    "TEVSACAS9": "50"
    }

class models:

    optimizer_selection=0
    nucleotide_dimensions=0
    compile_options={"optimizer": "unset", "loss": "unset"}
    
    def __init__(self, model_name, summary, nt_dims=4, optimizer="adam", learning_rate=0.0005, loss_function="mean_squared_error", drop_rate=0.3,CNN_filters=128, window_size=3, CNN_drop=0.3, conv1D_padding="same", CNN_dense1=128, CNN_dense2=64, maxpool1D_padding="same", RNN_size=128, RNN_dense1=128, RNN_dense2=64, CNN_RNN_drop=0.3):
        
        print("\nINITIALIZATION OF MODELS\n")
        try:
            models.optimizer_selection = models.get_optimizer(optimizer, learning_rate)
            models.compile_options={"optimizer": models.optimizer_selection, "loss": loss_function}
        except Exception as e:
            print(f"ERROR: Optimizer name {optimizer} not found or learning rate {learning_rate} not supported.\n{e}")
        try:
            self.nucleotide_dimensions=nt_dims
            self.model = self.build_crisprHAL_GPU(modelVersionInputLength[model_name][0], drop_rate, CNN_filters, window_size, CNN_drop, conv1D_padding, CNN_dense1, CNN_dense2, maxpool1D_padding, RNN_size, RNN_dense1, RNN_dense2, CNN_RNN_drop)
            if summary: self.model.summary()
            print("\n" + model_name + " model loaded.\n")
        except Exception as e:
            print(f"ERROR: Model name {model_name} not found.\n{e}")
    
    def get_optimizer(optimizer_name, learning_rate):
        if optimizer_name == 'adam':
            return k.optimizers.Adam(learning_rate=learning_rate, clipnorm=0.5)
        else:
            raise ValueError(f"Unsupported optimizer: {optimizer_name}, or unsupported learning rate: {learning_rate}")
    
    def build_crisprHAL_GPU(self, input_length, drop_rate=0.3, CNN_filters=128, window_size=3, CNN_drop=0.3, conv1D_padding="same", CNN_dense1=128, CNN_dense2=64, maxpool1D_padding="same", RNN_size=128, RNN_dense1=128, RNN_dense2=64, CNN_RNN_drop=0.3):
        i = k.Input(shape=(input_length,self.nucleotide_dimensions), name="Input")

        # Joint CNN Block
        c1 = Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c1")(i)
        l1 = LeakyReLU(name="l1")(c1)
        p1 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p1")(l1)
        dr1 = Dropout(CNN_drop, name="dr1")(p1)

        # Multi-layer CNN Branch
        c2 = Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c2")(dr1)
        l2 = LeakyReLU(name="l2")(c2)
        p2 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p2")(l2)
        dr2 = Dropout(CNN_drop, name="dr2")(p2)

        c3 = Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c3")(dr2)
        l3 = LeakyReLU(name="l3")(c3)
        p3 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p3")(l3)
        dr3 = Dropout(CNN_drop, name="dr3")(p3)
        
        c4 = Conv1D(CNN_filters,window_size,padding=conv1D_padding, name="c4")(dr3)
        l4 = LeakyReLU(name="l4")(c4)
        p4 = MaxPooling1D(pool_size=2,padding=maxpool1D_padding, name="p4")(l4)
        dr4 = Dropout(CNN_drop, name="dr4")(p4)

        f = Flatten(name="f")(dr4)
        
        d1 = Dense(CNN_dense1, name="d1")(f)
        ld1 = LeakyReLU(name="ld1")(d1)
        drd1 = Dropout(CNN_drop, name="drd1")(ld1)
        d2 = Dense(CNN_dense2, name="d2")(drd1)
        ld2 = LeakyReLU(name="ld2")(d2)
        drd2 = Dropout(CNN_drop, name="drd2")(ld2)
        x_1d_o = Dense(1, name="x_1d_o")(drd2)

        # CuDNN NVIDIA-GPU Optimized RNN: Bidirectional Gated Recurrent Unit Branch
        r1 = Bidirectional(GRU(RNN_size, kernel_initializer='he_normal', dropout=0.0, activation='tanh', recurrent_dropout=0.0, recurrent_activation='sigmoid', unroll=False, use_bias=True, reset_after=True), name="r1")(dr1)
        rd1 = Dropout(drop_rate, name="rd1")(r1)
        
        f_LSTM = Flatten(name="f_LSTM")(rd1)
        
        d1_LSTM = Dense(RNN_dense1, name="d1_LSTM")(f_LSTM)
        ld1_LSTM = LeakyReLU(name="ld1_LSTM")(d1_LSTM)
        drd1_LSTM = Dropout(CNN_RNN_drop, name="drd1_LSTM")(ld1_LSTM)
        d2_LSTM = Dense(RNN_dense2, name="d2_LSTM")(drd1_LSTM)
        ld2_LSTM = LeakyReLU(name="ld2_LSTM")(d2_LSTM)
        drd2_LSTM = Dropout(CNN_RNN_drop, name="drd2_LSTM")(ld2_LSTM)
        x_LSTM_o = Dense(1, name="x_LSTM_o")(drd2_LSTM)

        # Concatenated Output Layers
        o_c = concatenate([x_1d_o,x_LSTM_o], axis=1, name="o_c")
        o = Dense(1, activation="linear", name="o")(o_c)

        model = Model(inputs=i, outputs=o)
        model.compile(**models.compile_options)
        return model

    def train(self, x_train, y_train, epochs=10, batch_size=1024, verbose=1):
        try:
            self.model.fit(x_train, y_train, epochs=epochs, batch_size=batch_size, verbose=verbose)
        except Exception as e:
            print(f"ERROR: There was an issue during model training with X_train, y_train, or parameters: epochs={epochs}, batch_size={batch_size}, verbose={verbose}\n{e}")
    
    def predict(self, x, batch_size=2048, verbose=1):
        try:
            return self.model.predict(x, batch_size=batch_size, verbose=verbose)
        except Exception as e:
            print(f"ERROR: There was an issue during model prediction with the data provided or parameters: batch_size={batch_size}, verbose={verbose}\n{e}")

    def load_model(self, basename, path="."):
        path = path.strip("/")
        try:
            self.model = k.models.load_model(path + "/models/" + basename + ".keras")
            print(f"Loaded the model from the path: {path}/models/{basename}.keras")
        except Exception as e:
            print(f"ERROR: There was an issue loading the model: {path}/models/{basename}.keras\n{e}") 
    
    def compile(self):
        try:
            self.model.compile(**models.compile_options)
        except Exception as e:
            print(f"ERROR: Could not compile the model with compilation options: {str(self.compile_options)}\n{e}")
    
    def summary(self):
        self.model.summary()
    
    def clear_backend():
        k.backend.clear_session()
