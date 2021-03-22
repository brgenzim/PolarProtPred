from keras.models import load_model
from os import popen
import random
import lxml.etree as ET
from lxml import etree
import pickle
def fasta(inf):
	seqs={}
	s=""
	ID=""
	InFile=open(inf)
	while 1:
		line=InFile.readline()
		if line=="":
			seqs[ID]=s
			break
		if line[0]==">":
			if ID!="":
				seqs[ID]=s
				break		
			ID=line.strip()
			s=""
		else:
			s=s+line.replace(" ","").strip()
	InFile.close()
	return [s,ID]

def idr(inf):
	s=""
	InFile=open(inf)
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
	return s

def top(inf):
	s2=""	
	s=""
	doc=etree.parse(inf)
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
	return s

	

def lcr(inf):
	data2=""
	data=""
	a=open(inf)
	while 1:
		al=a.readline()
		if al=="":
			break
		if al[:12]=="{'proteins':":
			data=al.strip().replace("'",'"')
		elif al.find("sequence")!=-1:
			data2=al.strip().replace("'",'"')
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
	return "".join(s)

def model(loc):
	mods={}
	mods["model0A"]=load_model(loc+"/Models/model_0_A")
	mods["model1A"]=load_model(loc+"/Models/model_1_A")
	mods["model2A"]=load_model(loc+"/Models/model_2_A")
	mods["model3A"]=load_model(loc+"/Models/model_3_A")
	mods["model4A"]=load_model(loc+"/Models/model_4_A")
	mods["model0B"]=load_model(loc+"/Models/model_0_B")
	mods["model1B"]=load_model(loc+"/Models/model_1_B")
	mods["model2B"]=load_model(loc+"/Models/model_2_B")
	mods["model3B"]=load_model(loc+"/Models/model_3_B")
	mods["model4B"]=load_model(loc+"/Models/model_4_B")
	mods["model0P"]=load_model(loc+"/Models/model_0_P")
	mods["model1P"]=load_model(loc+"/Models/model_1_P")
	mods["model2P"]=load_model(loc+"/Models/model_2_P")
	mods["model3P"]=load_model(loc+"/Models/model_3_P")
	mods["model4P"]=load_model(loc+"/Models/model_4_P")
	mods["model0O"]=load_model(loc+"/Models/model_0_O")
	mods["model1O"]=load_model(loc+"/Models/model_1_O")
	mods["model2O"]=load_model(loc+"/Models/model_2_O")
	mods["model3O"]=load_model(loc+"/Models/model_3_O")
	mods["model4O"]=load_model(loc+"/Models/model_4_O")
	mods["model0Am"]=load_model(loc+"/Models/motifA_0")
	mods["model1Am"]=load_model(loc+"/Models/motifA_1")
	mods["model2Am"]=load_model(loc+"/Models/motifA_2")
	mods["model3Am"]=load_model(loc+"/Models/motifA_3")
	mods["model4Am"]=load_model(loc+"/Models/motifA_4")
	mods["model0Bm"]=load_model(loc+"/Models/motifB_0")
	mods["model1Bm"]=load_model(loc+"/Models/motifB_1")
	mods["model2Bm"]=load_model(loc+"/Models/motifB_2")
	mods["model3Bm"]=load_model(loc+"/Models/motifB_3")
	mods["model4Bm"]=load_model(loc+"/Models/motifB_4")
	mods["model0POm"]=load_model(loc+"/Models/motifPO_0")
	mods["model1POm"]=load_model(loc+"/Models/motifPO_1")
	mods["model2POm"]=load_model(loc+"/Models/motifPO_2")
	mods["model3POm"]=load_model(loc+"/Models/motifPO_3")
	mods["model4POm"]=load_model(loc+"/Models/motifPO_4")
	mods["final0"]=load_model(loc+"/Models/final0")
	mods["final1"]=load_model(loc+"/Models/final1")
	mods["final2"]=load_model(loc+"/Models/final2")
	mods["final3"]=load_model(loc+"/Models/final3")
	mods["final4"]=load_model(loc+"/Models/final4")
	mods["finalB0"]=load_model(loc+"/Models/finalB0")
	mods["finalB1"]=load_model(loc+"/Models/finalB1")
	mods["finalB2"]=load_model(loc+"/Models/finalB2")
	mods["finalB3"]=load_model(loc+"/Models/finalB3")
	mods["finalB4"]=load_model(loc+"/Models/finalB4")	
	return mods

def motif(loc):
	InFile = open(loc+"motif","rb")
	[mot,cat]=pickle.load(InFile)
	return [mot,cat]	
	
def species(loc):
	specs={}
	InFile=open(loc+"speindex.txt")
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
	return specs
def output_create(name):
	OutFile=open(name,"w")
	return OutFile

def output_write(OutFile,pred,score,name):
	OutFile.write(name[1:]+" "+pred+" "+str(score)+"\n")

	return ""

def output_close(OutFile):

	OutFile.close()
	return ""

def createtmpfile(loc,seq):
	Name=str(random.random())[2:]
	OutFile=open(loc+"/"+Name,"w")
	OutFile.write(">id\n"+seq)
	OutFile.close()
	return Name

def removetmpfile(name):
	tmp=popen("rm "+name)
	content=tmp.readlines()
	return ""




	


