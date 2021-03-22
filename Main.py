import ManageFiles
import Dependencies
import Tools
import sys
import hashlib

#Please edit these path
DIR=""
BLAST=""

#Nothing to do below
SeqFile=sys.argv[1]
TopologyFile=sys.argv[2]
DisorderedFile=sys.argv[3]
LowComplexityFile=sys.argv[4]
ResFile=sys.argv[5]
Mode=sys.argv[6]

print("Reading files")
[Sequence,Identifier]=ManageFiles.fasta(SeqFile)
Topology=ManageFiles.top(TopologyFile)
Disordered=ManageFiles.idr(DisorderedFile)
LCR=ManageFiles.lcr(LowComplexityFile)
Species=ManageFiles.species(DIR)
[Motifs,Sides]=ManageFiles.motif(DIR)
print("Reading model files")
Models=ManageFiles.model(DIR)
Output=ManageFiles.output_create(ResFile)

if Topology.find("M")==-1:
	print("Result")
	print("nonTM 1 "+key[1:])
	proceed=ManageFiles.output_write(Output,"nonTM","1",key,"","")		
else:
	print("Calculating sequence features")
	TmpFile=ManageFiles.createtmpfile(DIR,Sequence)
	Conservation=Dependencies.runblast(BLAST,DIR,TmpFile,Species,Sequence)
	print("Calculating featurematrix")
	Matrix=Tools.featuremx(Sequence,Topology,Disordered,LCR,Conservation)
	print("Scanning motifs")
	Vector=Tools.motif(Sequence,Topology,Disordered,Conservation,Motifs,Sides)
	print("Predicting localization")
	[Prediction,Score]=Tools.neural(Matrix,Vector,Models,Mode)
	print("Result")
	print(Identifier[1:]+" "+Prediction+" "+str(Score)+" ")
	proceed=ManageFiles.removetmpfile(TmpFile)
	proceed=ManageFiles.output_write(Output,Prediction,Score,Identifier)
print("Result is written in "+ResFile)
proceed=ManageFiles.output_close(Output)
print("Done")
