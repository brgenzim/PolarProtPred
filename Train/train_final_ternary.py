from __future__ import division 
from __future__ import print_function
import keras
import numpy as np
import re
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
from keras.models import load_model
from keras.callbacks import EarlyStopping
import lxml.etree as ET
from lxml import etree
import pickle

ABC="STPFILMVADERKHYGCNQW"

TRAINSEQUENCES="test.fas"#fasta format sequences
VALIDATIONSEQUENCES="test.fas"#fasta format sequences
SUBPREDICTORS="SUB"#result of subpredictors, where each line contains sequence indentifier, then first coordinate of various model outputs (order: CNN_Apical0-4,CNN_Basolateral0-4,CNN_Plasmal0-4,CNN_Other0-4,FCNN_Apical0-4,FCNN_Basolateral0-4,FCNN_Other0-4
LABELS="test.label"#each line contains sequence identifier and label (0/1) separated by tab/space


seqs={}
ID=""
s=""
inf=open(TRAINSEQUENCES)
while 1:
	l=inf.readline()
	if l=="":
		break
	if l[0]==">":
		ID=l[1:].strip()
	else:
		try:
			seqs[ID]=seqs[ID]+l.strip()
		except KeyError:
			seqs[ID]=l.strip()
inf.close()
labels={}
inf=open(LABELS)
while 1:
	l=inf.readline()
	if l=="":
		break
	l=l.split()
	labels[l[0]]=int(l[1].strip())
inf.close()

Y=[]
X=[]
for UP_ID in seqs:
	inf=open(SUBPREDICTORS+"/"+UP_ID+".result")
	while 1:
		l=inf.readline()
		if l=="":
			break
		l=l.split()
		X.append([])
		for i in range (1, len(l)):
			X[len(X)-1].append(float(l[i]))
		Y.append([0,0,0])
		Y[len(Y)-1][labels[l[0]]]=1


XT=np.array(X)
YT=np.array(Y)

seqs={}
ID=""
s=""
inf=open(VALIDATIONSEQUENCES)
while 1:
	l=inf.readline()
	if l=="":
		break
	if l[0]==">":
		ID=l[1:].strip()
	else:
		try:
			seqs[ID]=seqs[ID]+l.strip()
		except KeyError:
			seqs[ID]=l.strip()
inf.close()
Y2=[]
Y=[]
X=[]
for UP_ID in seqs:
	inf=open(SUBPREDICTORS+"/"+UP_ID+".result")
	while 1:
		l=inf.readline()
		if l=="":
			break
		l=l.split()
		X.append([])
		for i in range (1, len(l)):
			X[len(X)-1].append(float(l[i]))
		Y.append([0,0,0])
		Y[len(Y)-1][labels[l[0]]]=1
XV=np.array(X)
YV=np.array(Y)
modeldense = Sequential()			
modeldense.add(Dense(10,input_dim=35, activation='sigmoid'))
modeldense.add(Dense(3, activation='softmax'))
modeldense.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
modeldense.fit(XT, YT, batch_size=10, epochs=100, verbose=1, validation_data=(XV, YV))



