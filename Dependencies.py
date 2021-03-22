from os import popen
import lxml.etree as ET
from lxml import etree

def runblast(blast,loc,fasta,specs,seq):
	tmp=popen(blast+" -p blastp -d "+loc+"uniprot_sprot.fas -e 0.00001 -i "+fasta+" -m 7 -F F")
	content=tmp.read()
	doc = etree.fromstring(content)
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
	VECTOR.append(seq)
	for i in range (0, len(hitid)):
		if (hto[i]-hfrom[i])/len(seq)>0.75:
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
				for j in range (qto[i],len(seq)):
					s=s+"-"
				#print(s)
				VECTOR.append(s)

	string=""
	for i in range (0, len(seq)):
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
	return string
