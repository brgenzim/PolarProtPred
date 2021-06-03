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
import lxml.etree as ET
from lxml import etree
import pickle

ABC="STPFILMVADERKHYGCNQW"

TRAINSEQUENCES="test.fas"#fasta format sequences
VALIDATIONSEQUENCES="test.fas"#fasta format sequences
MOTIFS="motif"#dumped motifs, attached
SPECIES="speindex.txt"#speindex.txt from SwissProt
LCRDIR="LCR"#json formatted LCR files (standard Platoloco output)
TOPDIR="TOP"#xml formatted topology files (standard cctop output)
IDRDIR="IDR"#txt formatted disordered files (standard iupred output)
BLASTDIR="BLAST"#xml formatted blast files (standard blast -m 7 output)
LABELS="test.label"#each line contains sequence identifier and label (0/1) separated by tab/space

InFile = open(MOTIFS, 'rb')
[motifs,categ]=pickle.load(InFile)

specs={}
InFile=open(SPECIES)
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
line=InFile.readline()
while 1:
	line=InFile.readline()
	if line=="":
		break
	if line.find(")")!=-1:	
		if line.find(",")!=-1:
			line=line.replace("),",") ,")
			line=line.split()

			for i in range (0, len(line),3):
				try:
					specs[line[i+1].replace("(","").replace(")","")]=line[i].split("_")[1]
				except IndexError:
					pass
InFile.close()

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

lcrs={}
tops={}
idrs={}
cons={}
for UP_ID in seqs:
	inf=open(LCRDIR+"/"+UP_ID+".json")
	while 1:
		l=inf.readline()
		if l=="":
			break
		if l[:12]=="{'proteins':":
			data=l.strip().replace("'",'"')
		elif l.find("sequence")!=-1:
			data2=l.strip().replace("'",'"')
	s=[]
	proc_dat = eval(data) 
	proc_dat2 = eval(data2) 
	for key in proc_dat2:
		if key=="sequence":
			for i in range (0, len(proc_dat2[key])):
				s.append("0")
	for key in proc_dat:
		if key=="proteins":
			for key2 in proc_dat[key][0]:
				if key2=="SEG":
					for j in range (0, len(proc_dat[key][0][key2])):
						for k in range (int(proc_dat[key][0][key2][j][0])-1,int(proc_dat[key][0][key2][j][1])):
							s[k]="1"
	lcrs[UP_ID]="".join(s)
	inf.close()

	s2=""	
	s=""
	doc=etree.parse(TOPDIR+"/"+UP_ID+".xml")
	for elt2 in doc.getiterator():
		if elt2.tag=="Sequence":
			s=""
			for elt3 in elt2.getiterator():
				if elt3.tag=="Seq":
					s2=elt3.text.replace("\n","")

		if elt2.tag=="Topology":
			s=""
			for elt3 in elt2.getiterator():
				if elt3.tag=="Region":
				        for i in range (int(elt3.attrib["from"]),int(elt3.attrib["to"])+1):
				                s=s+elt3.attrib["loc"]
			
	if s=="":
		for i in range (0, len(s2)):
			s=s+"#"
	tops[UP_ID]=s


	s=""
	InFile=open(IDRDIR+"/"+UP_ID+".txt")
	while 1:
		line=InFile.readline()
		if line=="":
			break
		if line[0]!="#":
			line=line.split()
			if float(line[2])>=0.4:
				s=s+"1"
			else:
				s=s+"0"
	InFile.close()
	idrs[UP_ID]=s

	doc=etree.parse(BLASTDIR+"/"+UP_ID+".xml")
	hitid=[]
	qto=[]
	qfrom=[]
	hto=[]
	hfrom=[]
	midline=[]
	qseq=[]
	hseq=[]
	hitlen=[]
	ID=""
	hitlength=0
	for elt in doc.getiterator('Hit'):   
		bool1=False     
		bool2=False    
		for elt in elt.getiterator():
			if elt.tag=="Hit_def":
				ID=elt.text
				if ID!="id":
					hitid.append(ID)
					bool1=True
			if elt.tag=="Hit_len":
				if ID!="id":
					hitlength=int(elt.text)
					hitlen.append(int(elt.text))
					bool2=True
			if bool1==True and bool2==True:
				for elt in elt.getiterator('Hit_hsps'):
					for elt in elt.getiterator():
						if elt.tag=="Hsp_query-from":
							qfrom.append(int(elt.text))
							if len(qfrom)!=len(hitid):
								hitid.append(ID)  
							if len(qfrom)!=len(hitlen):	      
								hitlen.append(hitlength)
						if elt.tag=="Hsp_query-to":
							qto.append(int(elt.text))
						if elt.tag=="Hsp_hit-from":
							hfrom.append(int(elt.text))
						if elt.tag=="Hsp_hit-to":
							hto.append(int(elt.text))
						if elt.tag=="Hsp_qseq":
							qseq.append(elt.text)
						if elt.tag=="Hsp_hseq":
							hseq.append(elt.text)														   
						if elt.tag=="Hsp_midline":
							midline.append(elt.text)
	candidate=[]
	VECTOR=[]
	VECTOR.append(seqs[UP_ID])
	for i in range (0, len(hitid)):
		if (hto[i]-hfrom[i])/len(UP_ID)>0.75:
			bool1=True
			try:
			
				for j in range (0, len(candidate)):
					if candidate[j]==specs[hitid[i]]:
						bool1=False
			except KeyError:
				bool1=False
			if bool1==True:
				try:
					candidate.append(specs[hitid[i]])
				except KeyError:
					pass
				s=""
				for j in range (0, qfrom[i]-1):
					s=s+"-"
				for j in range (0, len(midline[i])):
					s=s+midline[i][j]	
				for j in range (qto[i],len(seqs[UP_ID])):
					s=s+"-"
				#print(s)
				VECTOR.append(s)

	string=""
	for i in range (0, len(seqs[UP_ID])):
		score=0
		hits=0
		for j in range (0, len(VECTOR)):
			if VECTOR[j][i]==" ":			

				hits+=1
			elif VECTOR[j][i]=="-":
				pass
			elif VECTOR[j][i]=="+":
				hits+=1
				score+=0.5
			else:
				hits+=1
				score+=1
		try:
			if str(score/hits)[0]=="1":
				string=string+"9"
			else:
				string=string+str(score/hits)[2]
		except Exception:
			string=string+"0"
	cons[UP_ID]=string




X=[]
Y=[]
for key in seqs:
	
	X.append([0,0,0,0,0,0,0])
	Y.append([])
	if labels[key]==1:
		Y[len(Y)-1].append(1)
		Y[len(Y)-1].append(0)

	else:
		Y[len(Y)-1].append(0)
		Y[len(Y)-1].append(1)
	for i in range (0, len(motifs)):	
		count=[0,0]#api in, api out, vasi in, basi out, other in, other out, other
		xx=seqs[key]
		p = re.compile(motifs[i])
		bool1=[False,False,False,False,False,False,False]
		for m in p.finditer(xx):
			b1=False

			if tops[key][m.start()]=="O" or tops[key][m.end()-1]=="O":
				if idrs[key][m.start():m.end()].find("1")!=-1:
					if cons[key][m.start():m.end()].find("6")!=-1 or cons[key][m.start():m.end()].find("7")!=-1 or cons[key][m.start():m.end()].find("8")!=-1 or cons[key][m.start():m.end()].find("9")!=-1:
						b1=True
	
			if b1==True:
				count[0]+=1
				if categ[motifs[i]].find("API")!=-1:
					X[len(X)-1][0]+=1
				if categ[motifs[i]].find("BAS")!=-1:
					X[len(X)-1][1]+=1
	
			if b1==False:
				if tops[key][m.start()]=="I" or tops[key][m.end()-1]=="I":
					if idrs[key][m.start():m.end()].find("1")!=-1:
						if cons[key][m.start():m.end()].find("6")!=-1 or cons[key][m.start():m.end()].find("7")!=-1 or cons[key][m.start():m.end()].find("8")!=-1 or cons[key][m.start():m.end()].find("9")!=-1:
							b1=True
	
				if b1==True:
					count[1]+=1
					if categ[motifs[i]].find("API")!=-1:
						X[len(X)-1][2]+=1
					if categ[motifs[i]].find("BAS")!=-1:
						X[len(X)-1][3]+=1

	countt=0
	counts=0
	count=0
	countp=0
	for i in range (0, len(seqs[key])):
		if tops[key][i]=="O":
			count+=1
			if seqs[key][i]=="S":
				counts+=1
			if seqs[key][i]=="T":
				countt+=1
			if seqs[key][i]=="P":
				countp+=1
	ins=0
	outs=0			
	for i in range (0, len(seqs[key])):
		if tops[key][i]=="O":
			outs+=1
		else:
			ins+=1


	if count>0:
		X[len(X)-1][4]=counts
		X[len(X)-1][5]=countt
		X[len(X)-1][6]=countp
	else:
		X[len(X)-1][4]=0
		X[len(X)-1][5]=0
		X[len(X)-1][6]=0

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
labels={}
inf=open(LABELS)
while 1:
	l=inf.readline()
	if l=="":
		break
	l=l.split()
	labels[l[0]]=int(l[1].strip())
inf.close()

lcrs={}
tops={}
idrs={}
cons={}
for UP_ID in seqs:
	inf=open(LCRDIR+"/"+UP_ID+".json")
	while 1:
		l=inf.readline()
		if l=="":
			break
		if l[:12]=="{'proteins':":
			data=l.strip().replace("'",'"')
		elif l.find("sequence")!=-1:
			data2=l.strip().replace("'",'"')
	s=[]
	proc_dat = eval(data) 
	proc_dat2 = eval(data2) 
	for key in proc_dat2:
		if key=="sequence":
			for i in range (0, len(proc_dat2[key])):
				s.append("0")
	for key in proc_dat:
		if key=="proteins":
			for key2 in proc_dat[key][0]:
				if key2=="SEG":
					for j in range (0, len(proc_dat[key][0][key2])):
						for k in range (int(proc_dat[key][0][key2][j][0])-1,int(proc_dat[key][0][key2][j][1])):
							s[k]="1"
	lcrs[UP_ID]="".join(s)
	inf.close()

	s2=""	
	s=""
	doc=etree.parse(TOPDIR+"/"+UP_ID+".xml")
	for elt2 in doc.getiterator():
		if elt2.tag=="Sequence":
			s=""
			for elt3 in elt2.getiterator():
				if elt3.tag=="Seq":
					s2=elt3.text.replace("\n","")

		if elt2.tag=="Topology":
			s=""
			for elt3 in elt2.getiterator():
				if elt3.tag=="Region":
				        for i in range (int(elt3.attrib["from"]),int(elt3.attrib["to"])+1):
				                s=s+elt3.attrib["loc"]
			
	if s=="":
		for i in range (0, len(s2)):
			s=s+"#"
	tops[UP_ID]=s


	s=""
	InFile=open(IDRDIR+"/"+UP_ID+".txt")
	while 1:
		line=InFile.readline()
		if line=="":
			break
		if line[0]!="#":
			line=line.split()
			if float(line[2])>=0.4:
				s=s+"1"
			else:
				s=s+"0"
	InFile.close()
	idrs[UP_ID]=s

	doc=etree.parse(BLASTDIR+"/"+UP_ID+".xml")
	hitid=[]
	qto=[]
	qfrom=[]
	hto=[]
	hfrom=[]
	midline=[]
	qseq=[]
	hseq=[]
	hitlen=[]
	ID=""
	hitlength=0
	for elt in doc.getiterator('Hit'):   
		bool1=False     
		bool2=False    
		for elt in elt.getiterator():
			if elt.tag=="Hit_def":
				ID=elt.text
				if ID!="id":
					hitid.append(ID)
					bool1=True
			if elt.tag=="Hit_len":
				if ID!="id":
					hitlength=int(elt.text)
					hitlen.append(int(elt.text))
					bool2=True
			if bool1==True and bool2==True:
				for elt in elt.getiterator('Hit_hsps'):
					for elt in elt.getiterator():
						if elt.tag=="Hsp_query-from":
							qfrom.append(int(elt.text))
							if len(qfrom)!=len(hitid):
								hitid.append(ID)  
							if len(qfrom)!=len(hitlen):	      
								hitlen.append(hitlength)
						if elt.tag=="Hsp_query-to":
							qto.append(int(elt.text))
						if elt.tag=="Hsp_hit-from":
							hfrom.append(int(elt.text))
						if elt.tag=="Hsp_hit-to":
							hto.append(int(elt.text))
						if elt.tag=="Hsp_qseq":
							qseq.append(elt.text)
						if elt.tag=="Hsp_hseq":
							hseq.append(elt.text)														   
						if elt.tag=="Hsp_midline":
							midline.append(elt.text)
	candidate=[]
	VECTOR=[]
	VECTOR.append(seqs[UP_ID])
	for i in range (0, len(hitid)):
		if (hto[i]-hfrom[i])/len(UP_ID)>0.75:
			bool1=True
			try:
			
				for j in range (0, len(candidate)):
					if candidate[j]==specs[hitid[i]]:
						bool1=False
			except KeyError:
				bool1=False
			if bool1==True:
				try:
					candidate.append(specs[hitid[i]])
				except KeyError:
					pass
				s=""
				for j in range (0, qfrom[i]-1):
					s=s+"-"
				for j in range (0, len(midline[i])):
					s=s+midline[i][j]	
				for j in range (qto[i],len(seqs[UP_ID])):
					s=s+"-"
				#print(s)
				VECTOR.append(s)

	string=""
	for i in range (0, len(seqs[UP_ID])):
		score=0
		hits=0
		for j in range (0, len(VECTOR)):
			if VECTOR[j][i]==" ":			

				hits+=1
			elif VECTOR[j][i]=="-":
				pass
			elif VECTOR[j][i]=="+":
				hits+=1
				score+=0.5
			else:
				hits+=1
				score+=1
		try:
			if str(score/hits)[0]=="1":
				string=string+"9"
			else:
				string=string+str(score/hits)[2]
		except Exception:
			string=string+"0"
	cons[UP_ID]=string




X=[]
Y=[]
for key in seqs:
	
	X.append([0,0,0,0,0,0,0])
	Y.append([])
	if labels[key]==1:
		Y[len(Y)-1].append(1)
		Y[len(Y)-1].append(0)

	else:
		Y[len(Y)-1].append(0)
		Y[len(Y)-1].append(1)
	for i in range (0, len(motifs)):	
		count=[0,0]#api in, api out, vasi in, basi out, other in, other out, other
		xx=seqs[key]
		p = re.compile(motifs[i])
		bool1=[False,False,False,False,False,False,False]
		for m in p.finditer(xx):
			b1=False

			if tops[key][m.start()]=="O" or tops[key][m.end()-1]=="O":
				if idrs[key][m.start():m.end()].find("1")!=-1:
					if cons[key][m.start():m.end()].find("6")!=-1 or cons[key][m.start():m.end()].find("7")!=-1 or cons[key][m.start():m.end()].find("8")!=-1 or cons[key][m.start():m.end()].find("9")!=-1:
						b1=True
	
			if b1==True:
				count[0]+=1
				if categ[motifs[i]].find("API")!=-1:
					X[len(X)-1][0]+=1
				if categ[motifs[i]].find("BAS")!=-1:
					X[len(X)-1][1]+=1
	
			if b1==False:
				if tops[key][m.start()]=="I" or tops[key][m.end()-1]=="I":
					if idrs[key][m.start():m.end()].find("1")!=-1:
						if cons[key][m.start():m.end()].find("6")!=-1 or cons[key][m.start():m.end()].find("7")!=-1 or cons[key][m.start():m.end()].find("8")!=-1 or cons[key][m.start():m.end()].find("9")!=-1:
							b1=True
	
				if b1==True:
					count[1]+=1
					if categ[motifs[i]].find("API")!=-1:
						X[len(X)-1][2]+=1
					if categ[motifs[i]].find("BAS")!=-1:
						X[len(X)-1][3]+=1

	countt=0
	counts=0
	count=0
	countp=0
	for i in range (0, len(seqs[key])):
		if tops[key][i]=="O":
			count+=1
			if seqs[key][i]=="S":
				counts+=1
			if seqs[key][i]=="T":
				countt+=1
			if seqs[key][i]=="P":
				countp+=1
	ins=0
	outs=0			
	for i in range (0, len(seqs[key])):
		if tops[key][i]=="O":
			outs+=1
		else:
			ins+=1


	if count>0:
		X[len(X)-1][4]=counts
		X[len(X)-1][5]=countt
		X[len(X)-1][6]=countp
	else:
		X[len(X)-1][4]=0
		X[len(X)-1][5]=0
		X[len(X)-1][6]=0

XV=np.array(X)
YV=np.array(Y)

modeldense = Sequential()			
modeldense.add(Dense(28,input_dim=7, activation='sigmoid'))
modeldense.add(Dropout(0.5))
modeldense.add(Dense(2, activation='softmax'))
modeldense.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
modeldense.fit(XT, YT, batch_size=10, epochs=100, verbose=1,validation_data=(XV, YV))





