import re
import getopt
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import Bio
from Bio import SeqIO
def usage():
	print ("--in *.Units.fa")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	name=''
	start=0
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		elif opt in ("--name"):
			name=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,name

def uniqUnits(pwd1,pwd2,name):
	#pwd1,name=getoptions()

	pwdo=name+'.Units.uniq.fa'
	fileo=open(pwdo,'w')
	
	dict0={}
	for gseq in SeqIO.parse(pwd2,'fasta'):
		dict0[gseq.id]=len(gseq.seq)


	dict1={}

	for gseq in SeqIO.parse(pwd1,'fasta'):
		ID=gseq.id.split('.unit-')[0]
		len1=dict0[ID]
		if len(gseq.seq)/float(len1)>=1.5:
			continue
		if ID in dict1:
			list1=dict1[ID]
		else:
			list1=[]
		len1=gseq.description.split('\t')[-1].split('-')
		s=int(len1[0])
		e=int(len1[1])
		list1.append([s,e,gseq])
		dict1[ID]=list1
	lenall=0
	for k,v in dict1.items():
		list1=sorted(v,key=(lambda x:x[0]))
		l1=list1[0]
		s,e,Seq=l1
		for l2 in list1[1:]:
			s,e,Seq=l1
			s1,e1,Seq1=l2
			sl1=len(Seq.seq)
			sl2=len(Seq1.seq)
			x=range(s-1,e,1)
			y=range(s1-1,e1,1)
			t=set(x).intersection(set(y))
			b=min([sl1,sl2])
			a=float(len(t))/float(b)
			if a>0.5:
				if sl2<sl1:
					l1=l2
			else:
				l='>'+Seq.description+"\n"+Seq.seq+"\n"
				fileo.writelines(l)
				lenall+=len(Seq.seq)
				l1=l2
		s,e,Seq=l1
		l='>'+Seq.description+"\n"+Seq.seq+"\n"
		fileo.writelines(l)
		lenall+=len(Seq.seq)
	print (lenall)
	print ('welldone')
	fileo.close()
	pwdo1=name+'.Consensus.uniq.fa'
	fileo=open(pwdo1,'w')
	IDlist=[]
	for gseq in SeqIO.parse(pwdo,'fasta'):
		ID=gseq.id.split('.unit-')[0]
		if ID not in IDlist:
			IDlist.append(ID)
	for gseq in SeqIO.parse(pwd2,'fasta'):
		if gseq.id in IDlist:
			l='>'+gseq.description+"\n"+gseq.seq+"\n"
			fileo.writelines(l)
	fileo.close()
