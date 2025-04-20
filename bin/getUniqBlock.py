import re
import Bio
from Bio import SeqIO

import getopt
import sys
def usage():
	print ("--in *.Consensus.fa")
	print ("--unit *.Units.fa")
	print ("--stat *.stat.txt")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","unit=","stat="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	stat=''
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
		if opt in ("--stat"):
			stat=value
		elif opt in ("--name"):
			name=value	
		elif opt in ("-h"):
			usage()
			sys.exit()
	getUniqBlock(seq1,seq2,stat,name)
#	return seq1,seq2,stat,name
def getUniqBlock(pwd1,pwd2,pwd3,name):
	#pwd1,pwd2,pwd3,name=getoptions()
	pwdo=name+'.block.fa'
	pwdo1=name+'.block.units.fa'
	fileo=open(pwdo,'w')
	fileo1=open(pwdo1,'w')
	
	file1=open(pwd3,'r')
	dict1={}
	for row in file1:
		row1=row.rstrip().split('\t')
		u=row1[-1].split(';')
		dict1[row1[0]]=[int(row1[1]),u]
	file1.close()
	print (len(dict1))
	dict2={}
	for gseq in SeqIO.parse(pwd2,'fasta'):
		ID=gseq.id.split('.unit-')[0]
		loc=gseq.description.split('\t')[-1].split('-')
		if ID in dict2:
			list1=dict2[ID]
		else:
			list1=[]
		list1.append(int(loc[0]))
		list1.append(int(loc[1]))
		dict2[ID]=list1
	list11=[]
	print (len(dict2))
	for gseq in SeqIO.parse(pwd1,'fasta'):
		if gseq.id not in dict2:
			continue
		i=gseq.description.split('\t')
		s=min(dict2[gseq.id])
		e=max(dict2[gseq.id])
		cp=len(dict1[gseq.id][1])
		list11.append([gseq.id,s,e,cp,gseq])
	
	list1=sorted(list11,key=(lambda x:x[1]))
	
	l1=list1[0]
	ID,s,e,cp,Seq=l1
	lenall=0
	IDlistunit=[]
	for l2 in list1[1:]:
		ID,s,e,cp,Seq=l1
		ID1,s1,e1,cp1,Seq1=l2
		x=range(s,e+1,1)
		y=range(s1,e1+1,1)
		t=set(x).intersection(set(y))
		if len(t)!=0:
			if cp1>=cp:
			#if (e1-s1)>=(e-s):
				l1=l2
		else:
			Seqo=Seq.description.split('\t')
			Seqo[2]=str(s)
			Seqo[3]=str(e)
			cp=str(len(Seq.seq))+"-"+str(len(dict1[Seq.id][1]))
			Seqo[1]=cp
			Seqo1="\t".join(Seqo)
			l='>'+Seqo1+"\n"+Seq.seq+"\n"
			fileo.writelines(l)
			v=dict1[Seq.id][1]
			IDlistunit=IDlistunit+v
			lenall+=len(Seq.seq)
			l1=l2
	ID,s,e,cp,Seq=l1
	l='>'+Seq.description+"\n"+Seq.seq+"\n"
	v=dict1[Seq.id][1]
	IDlistunit=IDlistunit+v
	fileo.writelines(l)
	lenall+=len(Seq.seq)
	fileo.close()
	j=0
	for gseq in SeqIO.parse(pwd2,'fasta'):
		if gseq.id in IDlistunit:
			l='>'+gseq.description+"\n"+gseq.seq+"\n"
			fileo1.writelines(l)
			j+=1
	fileo1.close()
	if len(IDlistunit)!=j:
		print ('error')
		print (IDlistunit[0])
		print (len(IDlistunit),j)
		print ('error')
	print ('welldone')
#getoptions()
