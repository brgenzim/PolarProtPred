import keras
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
import re

def featuremx(s,t,i,x,c):
	ABC="STPFILMVADERKHYGCNQW"
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
	return mx

def motif(s,t,i,c,mot,cat):
	X=[]
	X.append([0,0,0,0,0,0,0])
	for k in range (0, len(mot)):	
		p = re.compile(mot[k])
		for m in p.finditer(s):
			b1=False
			if t[m.start()]=="O" or t[m.end()-1]=="O":
				if i[m.start():m.end()].find("1")!=-1:
					if c[m.start():m.end()].find("6")!=-1:				
						b1=True
					elif c[m.start():m.end()].find("7")!=-1:
						b1=True
					elif c[m.start():m.end()].find("8")!=-1:
						b1=True
					elif c[m.start():m.end()].find("9")!=-1:
						b1=True
			if b1==True:
				if cat[mot[k]].find("API")!=-1:
					X[len(X)-1][0]+=1
				if cat[mot[k]].find("BAS")!=-1:
					X[len(X)-1][1]+=1
			if b1==False:
				if t[m.start()]=="I" or t[m.end()-1]=="I":
					if i[m.start():m.end()].find("1")!=-1:
						if c[m.start():m.end()].find("6")!=-1:				
							b1=True
						elif c[m.start():m.end()].find("7")!=-1:
							b1=True
						elif c[m.start():m.end()].find("8")!=-1:
							b1=True
						elif c[m.start():m.end()].find("9")!=-1:
							b1=True
				if b1==True:
					if cat[mot[k]].find("API")!=-1:
						X[len(X)-1][2]+=1
					if cat[mot[k]].find("BAS")!=-1:
						X[len(X)-1][3]+=1

	countt=0
	counts=0
	countp=0
	for k in range (0, len(s)):
		if t[k]=="O":
			if s[k]=="S":
				counts+=1
			if s[k]=="T":
				countt+=1
			if s[k]=="P":
				countp+=1
	X[len(X)-1][4]=counts
	X[len(X)-1][5]=countt
	X[len(X)-1][6]=countp
	X=np.array(X)
	return X
def neural(mx,vec,mod,pr):
	mat=[]
	mat.append(mx)
	mat=np.array(mat)
	mat = mat.reshape(mat.shape[0], 20, 20, 2)
	mat = mat.astype('float32')
	mat /= 255
	X=[]
	X.append(mod["model0Am"].predict(vec)[0][0])
	X.append(mod["model1Am"].predict(vec)[0][0])
	X.append(mod["model2Am"].predict(vec)[0][0])
	X.append(mod["model3Am"].predict(vec)[0][0])
	X.append(mod["model4Am"].predict(vec)[0][0])
	X.append(mod["model0Bm"].predict(vec)[0][0])
	X.append(mod["model1Bm"].predict(vec)[0][0])
	X.append(mod["model2Bm"].predict(vec)[0][0])
	X.append(mod["model3Bm"].predict(vec)[0][0])
	X.append(mod["model4Bm"].predict(vec)[0][0])
	X.append(mod["model0POm"].predict(vec)[0][0])
	X.append(mod["model1POm"].predict(vec)[0][0])
	X.append(mod["model2POm"].predict(vec)[0][0])
	X.append(mod["model3POm"].predict(vec)[0][0])
	X.append(mod["model4POm"].predict(vec)[0][0])
	X.append(mod["model0A"].predict(mat)[0][0])
	X.append(mod["model1A"].predict(mat)[0][0])
	X.append(mod["model2A"].predict(mat)[0][0])
	X.append(mod["model3A"].predict(mat)[0][0])
	X.append(mod["model4A"].predict(mat)[0][0])
	X.append(mod["model0B"].predict(mat)[0][0])
	X.append(mod["model1B"].predict(mat)[0][0])
	X.append(mod["model2B"].predict(mat)[0][0])
	X.append(mod["model3B"].predict(mat)[0][0])
	X.append(mod["model4B"].predict(mat)[0][0])
	X.append(mod["model0P"].predict(mat)[0][0])
	X.append(mod["model1P"].predict(mat)[0][0])
	X.append(mod["model2P"].predict(mat)[0][0])
	X.append(mod["model3P"].predict(mat)[0][0])
	X.append(mod["model4P"].predict(mat)[0][0])
	X.append(mod["model0O"].predict(mat)[0][0])
	X.append(mod["model1O"].predict(mat)[0][0])
	X.append(mod["model2O"].predict(mat)[0][0])
	X.append(mod["model3O"].predict(mat)[0][0])
	X.append(mod["model4O"].predict(mat)[0][0])
	X2=[]
	X2.append(X)
	X2=np.array(X2)
	if pr=="binary":
		xx0=mod["finalB0"].predict(X2)
		xx1=mod["finalB1"].predict(X2)
		xx2=mod["finalB2"].predict(X2)
		xx3=mod["finalB3"].predict(X2)
		xx4=mod["finalB4"].predict(X2)
		a=xx0[0][0]+xx1[0][0]+xx2[0][0]+xx3[0][0]+xx4[0][0]
		bop=xx0[0][1]+xx1[0][1]+xx2[0][1]+xx3[0][1]+xx4[0][1]
		if a>bop:
			return ["Apical",str(a/5)]
		else:
			return ["Basolateral",str(bop/5)]
	else:
		xx0=mod["final0"].predict(X2)
		xx1=mod["final1"].predict(X2)
		xx2=mod["final2"].predict(X2)
		xx3=mod["final3"].predict(X2)
		xx4=mod["final4"].predict(X2)

		a=xx0[0][0]+xx1[0][0]+xx2[0][0]+xx3[0][0]+xx4[0][0]
		b=xx0[0][1]+xx1[0][1]+xx2[0][1]+xx3[0][1]+xx4[0][1]
		op=xx0[0][2]+xx1[0][2]+xx2[0][2]+xx3[0][2]+xx4[0][2]
		if a==max(a,b,op):
			return ["Apical",str(a/5)]
		elif b==max(a,b,op):
			return ["Basolateral",str(b/5)]
		else:
			return ["Other",str(op/5)]

