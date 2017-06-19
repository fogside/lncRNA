import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from copy import deepcopy

from keras.models import Sequential
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import GlobalAveragePooling1D
from keras.layers.pooling import MaxPooling1D
from keras.layers import BatchNormalization
from keras.preprocessing import sequence
from keras.layers import Embedding
from keras.layers import Dense


MAX_SEQUENCE_LENGTH = 5000
N = 4000 #  количество сэмплов/2 в обучающей и валидирующей выборке
EMBEDDING_DIM = 100 # это число нельзя менять

#### Read fasta ####
def read_FASTA(file_name):
    with open(file_name, "r") as fn:
        text = fn.read().split(">")
    text = [x.split("\n") for x in text if x != ""]
    text = [[x[0],"".join(x[1:]).upper()] for x in text]
    text_dict = {line[0].split('|')[0]:line[1] for line in text}
    return text_dict


code = read_FASTA("../../fasta_data/train_fasta_gencode.fasta")
noncode = read_FASTA("../../fasta_data/train_fasta_lncPedia.fasta")
embedding_matrix = pd.read_csv("../data/embedding_matrix.csv", index_col=0)
encoder_dict = dict(zip(embedding_matrix.index, embedding_matrix.code))
embedding_matrix.index = embedding_matrix.code
embedding_matrix.drop('code', axis=1, inplace=True)

def make_kmers(seq, kmer=8, max_length=7000):
    res = []
    i=0
    while (i<(len(seq)-kmer)) and (i<max_length):
        tmp = seq[i:i+kmer]
        if 'N' in tmp:
            i+=1
            continue
        res.append(encoder_dict[tmp])
        i+=1
    return res

y = np.array([1]*N + [0]*N)
X_tmp = deepcopy(list(code.items())[:N])
X_tmp.extend(list(noncode.items())[:N])

X = np.array([make_kmers(v, max_length=MAX_SEQUENCE_LENGTH) for k,v in X_tmp])
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

embedding_layer = Embedding(len(embedding_matrix),
                            EMBEDDING_DIM,
                            weights=[np.array(embedding_matrix)],
                            input_length=MAX_SEQUENCE_LENGTH,
                            trainable=False, )

print("train: ", X_train.shape)
print("test: ", X_test.shape)
X_train = sequence.pad_sequences(X_train, maxlen=MAX_SEQUENCE_LENGTH)
X_test = sequence.pad_sequences(X_test, maxlen=MAX_SEQUENCE_LENGTH)

model = Sequential()
model.add(embedding_layer)
model.add(BatchNormalization())
model.add(Conv1D(filters=400, kernel_size=5, padding='same', activation='relu'))
model.add(GlobalAveragePooling1D())
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['acc'])
print(model.summary())
model.fit(X_train, y_train, epochs=15, batch_size=80, validation_data=(X_test, y_test))

scores = model.evaluate(X_test, y_test, verbose=0)
print("Accuracy: %.2f%%" % (scores[1]*100))


model.save("model_conv1D_400_v2.h5")
