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
SPECIES="speindex.txt"#speindex.txt from SwissProt
LCRDIR="LCR"#json formatted LCR files (standard Platoloco output)
TOPDIR="TOP"#xml formatted topology files (standard cctop output)
IDRDIR="IDR"#txt formatted disordered files (standard iupred output)
BLASTDIR="BLAST"#xml formatted blast files (standard blast -m 7 output)
LABELS="test.label"#each line contains sequence identifier and label (0/1) separated by tab/space

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
	
	s=seqs[key]
	t=tops[key]
	i=idrs[key]
	x=lcrs[key]
	c=cons[key]

	
	distance=[]
	for k in range (0, len(s)):
		if t[k]!="M":
			pos=k
			d1=-1
			d2=-1
			while 1:
				pos+=1
				if pos>len(t)-1:
					break
				if t[pos]=="M":
					d1=abs(pos-k)
					break
			pos=k

			while 1:
				pos-=1
				if pos<1:
					break
				if t[pos]=="M":
					d2=abs(pos-k)
					break
			if d1==-1:
				distance.append(d2)
			elif d2==-1:
				distance.append(d1)
			else:
				distance.append(min(d1,d2))
		else:
			distance.append(0)

	mx=[]
	for k in range (0, 20):
		mx.append([])
		for l in range (0, 20):
			mx[k].append([])
			for m in range (0,2):
				mx[k][l].append(0)

	for k in range (0, len(s)-1):
		if i[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][0]+=2*float(c[k])
				else:
					mx[p2][p1][0]+=2*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][0]+=2*float(c[k])
				else:
					mx[p2][p1][0]+=2*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+1]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][0]+=2*float(c[k])
			else:
				mx[p2][p1][0]+=2*float(c[k])
	for k in range (0, len(s)-2):
		if i[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][0]+=1*float(c[k])
				else:
					mx[p2][p1][0]+=1*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][0]+=1*float(c[k])
				else:
					mx[p2][p1][0]+=1*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+2]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][0]+=1*float(c[k])
			else:
				mx[p2][p1][0]+=1*float(c[k])
	for k in range (0, len(s)-1):
		if x[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][1]+=2*float(c[k])
				else:
					mx[p2][p1][1]+=2*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][1]+=2*float(c[k])
				else:
					mx[p2][p1][1]+=2*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+1]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][1]+=2*float(c[k])
			else:
				mx[p2][p1][1]+=2*float(c[k])
	for k in range (0, len(s)-2):
		if x[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][1]+=1*float(c[k])
				else:
					mx[p2][p1][1]+=1*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][1]+=1*float(c[k])
				else:
					mx[p2][p1][1]+=1*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+2]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][1]+=1*float(c[k])
			else:
				mx[p2][p1][1]+=1*float(c[k])
	for k in range (0, len(mx)):
		for l in range (0, len(mx[k])):
			if mx[k][l][0]>255:
				mx[k][l][0]=255
			if mx[k][l][1]>255:
				mx[k][l][1]=255
	X.append(mx)
	Y.append(labels[key])


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
	
	s=seqs[key]
	t=tops[key]
	i=idrs[key]
	x=lcrs[key]
	c=cons[key]

	
	distance=[]
	for k in range (0, len(s)):
		if t[k]!="M":
			pos=k
			d1=-1
			d2=-1
			while 1:
				pos+=1
				if pos>len(t)-1:
					break
				if t[pos]=="M":
					d1=abs(pos-k)
					break
			pos=k

			while 1:
				pos-=1
				if pos<1:
					break
				if t[pos]=="M":
					d2=abs(pos-k)
					break
			if d1==-1:
				distance.append(d2)
			elif d2==-1:
				distance.append(d1)
			else:
				distance.append(min(d1,d2))
		else:
			distance.append(0)

	mx=[]
	for k in range (0, 20):
		mx.append([])
		for l in range (0, 20):
			mx[k].append([])
			for m in range (0,2):
				mx[k][l].append(0)

	for k in range (0, len(s)-1):
		if i[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][0]+=2*float(c[k])
				else:
					mx[p2][p1][0]+=2*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][0]+=2*float(c[k])
				else:
					mx[p2][p1][0]+=2*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+1]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][0]+=2*float(c[k])
			else:
				mx[p2][p1][0]+=2*float(c[k])
	for k in range (0, len(s)-2):
		if i[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][0]+=1*float(c[k])
				else:
					mx[p2][p1][0]+=1*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][0]+=1*float(c[k])
				else:
					mx[p2][p1][0]+=1*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+2]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][0]+=1*float(c[k])
			else:
				mx[p2][p1][0]+=1*float(c[k])
	for k in range (0, len(s)-1):
		if x[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][1]+=2*float(c[k])
				else:
					mx[p2][p1][1]+=2*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+1]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][1]+=2*float(c[k])
				else:
					mx[p2][p1][1]+=2*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+1]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][1]+=2*float(c[k])
			else:
				mx[p2][p1][1]+=2*float(c[k])
	for k in range (0, len(s)-2):
		if x[k]=="1" and distance[k]>0:
			if t[k]=="O":					
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1>=p2:
					mx[p1][p2][1]+=1*float(c[k])
				else:
					mx[p2][p1][1]+=1*float(c[k])
			else:
				p1=0
				p2=0
				for l in range (0, len(ABC)):
					if s[k]==ABC[l]:
						p1+=l
					if s[k+2]==ABC[l]:
						p2+=l
				if p1<p2:
					mx[p1][p2][1]+=1*float(c[k])
				else:
					mx[p2][p1][1]+=1*float(c[k])
		else:
			p1=0
			p2=0
			for l in range (0, len(ABC)):
				if s[k]==ABC[l]:
					p1+=l
				if s[k+2]==ABC[l]:
					p2+=l
			if p1<p2:
				mx[p1][p2][1]+=1*float(c[k])
			else:
				mx[p2][p1][1]+=1*float(c[k])
	for k in range (0, len(mx)):
		for l in range (0, len(mx[k])):
			if mx[k][l][0]>255:
				mx[k][l][0]=255
			if mx[k][l][1]>255:
				mx[k][l][1]=255
	X.append(mx)
	Y.append(labels[key])

XV=np.array(X)
YV=np.array(Y)



num_classes=2
YT = keras.utils.to_categorical(YT, num_classes)
YV =  keras.utils.to_categorical(YV, num_classes)
input_shape = (20, 20, 2)
model = Sequential()
model.add(Conv2D(32, kernel_size=(3, 3),
		 activation='relu',
		 input_shape=input_shape))
model.add(Conv2D(64, (3, 3), activation='relu'))
model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='softmax'))
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1,patience=3)
model.compile(loss=keras.losses.categorical_crossentropy,optimizer=keras.optimizers.Adadelta(), metrics=['accuracy'])
model.fit(XT, YT, batch_size=10, epochs=20, verbose=1, validation_data=(XV, YV))






