import re
import Bio
from Bio import SeqIO

import getopt
import sys
def usage():
	print ("--in *.Consensus.fa")
	print ("--unit *.Units.fa")
	print ("--type *.stat.txt")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","unit=","type="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	type1=''
	name=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		if opt in ("--unit"):
			seq2=value
		if opt in ("--type"):
			type1=value
		elif opt in ("--name"):
			name=value	
		elif opt in ("-h"):
			usage()
			sys.exit()
	getTypeFromMerge(seq1,seq2,type1,name)
	return seq1,seq2,type1,name
def getTypeFromMerge(pwd1,pwd2,pwd3,name):
	pwdo=name+'.Type.fa'
	pwdo1=name+'.Type.units.fa'
	fileo=open(pwdo,'w')
	fileo1=open(pwdo1,'w')
	
	file1=open(pwd3,'r')
	dict1={}
	dictt={}
	for row in file1:
		row1=row.rstrip().split('\t')
		u=row1[-1].split(';')
		dict1[row1[0]]=row1[1]
		for i in u:
			dictt[i]=row1[0]
	file1.close()
	print (len(dict1))
	dict2={}
	for gseq in SeqIO.parse(pwd2,'fasta'):
		ID=gseq.id.split('.unit-')[0]
		if ID not in dictt:
			continue
		ty=dictt[ID]
		if ty in dict2:
			listt=dict2[ty]
		else:
			listt=[]
		listt.append(gseq)
		dict2[ty]=listt
	dict3={}
	dict31={}
	for gseq in SeqIO.parse(pwd1,'fasta'):
		if gseq.id not in dictt:
			continue
		ty=dictt[gseq.id]
		if ty in dict3:
			listt=dict3[ty]
		else:
			listt=[]
		listt.append(gseq)
		dict3[ty]=listt
		dict31[gseq.id]=gseq
	for k,v in dict1.items():
		Units=dict2[k]
		Cons1=dict1[k]
		if len(Units)>10:
			for Seq in Units:
				v=Seq.description.split('\t')
				l='>'+Seq.id+"\t"+str(k)+"\t"+"\t".join(v[1:])+"\n"+Seq.seq+"\n"
				fileo1.writelines(l)
			Cons=dict31[Cons1]
			v=Cons.description.split('\t')
			l='>'+Cons.id+"\t"+str(len(Units))+"\t"+"\t".join(v[1:])+"\t"+str(k)+"\n"+Cons.seq+"\n"
			fileo.writelines(l)
	fileo.close()
	fileo1.close()
	print ('welldone')
#getoptions()
