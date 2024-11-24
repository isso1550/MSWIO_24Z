# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 07:18:08 2022

@author: wykladowca
"""

#biblioteki i funkcje
SEED_NUM=0
import os 
os.environ['PYTHONHASHSEED']=str(SEED_NUM)


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
#

from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation,advanced_activations,LSTM
from keras.activations import relu, elu, sigmoid
from keras.optimizers import Adam, Nadam, RMSprop, Adagrad, Adadelta
from keras.callbacks import ModelCheckpoint,EarlyStopping
import tensorflow as tf
import random

def write(frame,name):
    writer=pd.ExcelWriter(name+'.xlsx')
    frame.to_excel(writer)
    writer.save()


def reset_random_seeds():
    os.environ['PYTHONHASHSEED']=str(SEED_NUM)
    tf.random.set_seed(SEED_NUM)    
    np.random.seed(SEED_NUM)
    random.seed(SEED_NUM)


###############################################################
#wczytywanie przykładowej ramki danych
data=pd.read_excel('df_ex.xlsx')

#Nomalizacja zmiennych
data['A']=data['A']/100
data['B']=data['B']/100
data['Res']=data['Res']/10000

#podział danych na treningowe/walidaycjne/testowe
df_tr_in=data.loc[0:74,['A','B']]
df_val_in=data.loc[75:100,['A','B']]
df_test_in=data.loc[101:200,['A','B']]

df_tr_out=data.loc[0:74,'Res']
df_val_out=data.loc[75:100,'Res']
df_test_out=data.loc[101:200,'Res']

#przeformatowanie danych na postać macierzy 3D do lstm
tr_samples_num=len(df_tr_in)
in_tr=df_tr_in.values.reshape(tr_samples_num,1,len(df_tr_in.columns)) #najpierw liczba próbek, potem liczba wierszy, potem kolumn
out_tr=df_tr_out.values.reshape(tr_samples_num,1,1)

val_samples_num=len(df_val_in)
in_val=df_val_in.values.reshape(val_samples_num,1,len(df_val_in.columns)) #najpierw liczba próbek, potem liczba wierszy, potem kolumn
out_val=df_val_out.values.reshape(val_samples_num,1,1)

test_samples_num=len(df_test_in)
in_test=df_test_in.values.reshape(test_samples_num,1,len(df_test_in.columns)) #najpierw liczba próbek, potem liczba wierszy, potem kolumn
out_test=df_test_out.values.reshape(test_samples_num,1,1)


#w razie gdyby tensorflow nie chciał przyjąć  macierzy 3d z numpy konwersja na natywne tensory
#in_tr=tf.convert_to_tensor(in_tr)
#out_tr=tf.convert_to_tensor(out_tr)



reset_random_seeds()

model = Sequential()
#(timesteps, data_dim)
model.add(LSTM(8,return_sequences=True, input_shape=(in_tr.shape[1],in_tr.shape[2]),activation='tanh'))
#model.add(Dropout(0.8))

#model.add(Dense(10, activation='relu'))
model.add(LSTM(4,return_sequences=True,activation='tanh'))
#model.add(LSTM(3,return_sequences=True,activation='sigmoid'))
#model.add(Dropout(0.2))
model.add(Dense(1, activation='linear'))

#opt=  Adam(lr=1e-4, decay=1e-5)
#sgd = SGD(lr=0.001, clipvalue=0.5,nesterov=True,momentum=0.8,clipnorm=1.0)
#opt = SGD(lr=0.001, decay=1e-6, momentum=0.8, nesterov=True)
#opt = Adagrad(lr=0.01, epsilon=1e-6, decay=0.0) #ładny wykres loss, mae ~0.115
#opt = RMSprop(lr=0.001, rho=0.9, epsilon=1e-6)
#opt =Adadelta(lr=1.0, rho=0.95, epsilon=None, decay=0.0)
opt= Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-07, decay=0.0001, amsgrad=True)
#opt= Nadam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-07)
#opt=Ftrl(lr=0.001,learning_rate_power=-0.5,initial_accumulator_value=0.1,l1_regularization_strength=0.0,l2_regularization_strength=0.0,l2_shrinkage_regularization_strength=0.0,)
model.compile(loss='mae',optimizer=opt)

pat=20
epo=200

early_stopping =EarlyStopping(monitor='val_loss', patience=pat)
#instrukja uczenia modelu
history=model.fit(in_tr, out_tr,validation_data=(in_val,out_val),epochs=epo,batch_size=128,verbose=1,shuffle=True,callbacks=[early_stopping]) 
#instrukja walidacji modelu
score = model.evaluate(in_test, out_test, batch_size=128) #instrukcja testowania modelu

ucz=history.history['loss'][-1]
spr=history.history['val_loss'][-1]
prognosis=model.predict(in_test, batch_size=None, verbose=0, steps=None) #wpisać jak chce się wartości przewidywane uzyskać
print(f'MAE test [%] = {score*100} ,uczenie= {ucz*100}, sprawdz={spr*100}')


prognosis2=prognosis.flatten().tolist()
out_test2=out_test.flatten().tolist()
out=pd.DataFrame(list(zip(prognosis2,out_test2)))
out.columns=['prognosis','reals']

write(out,'prognoza_iter')

plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Val'], loc='upper left')
plt.show()
