import re
import Bio
from Bio import SeqIO

import getopt
import sys
def usage():
	print ("--in *.Consensus.fa")
	print ("--unit *.Units.fa")
	print ("--TE *.TE.out")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","unit=","TE=","start=","end="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	TE=''
	name=''
	start=''
	end=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		if opt in ("--unit"):
			seq2=value
		if opt in ("--TE"):
			TE=value
		elif opt in ("--name"):
			name=value
		elif opt in ("--start"):
			start=int(value)
		elif opt in ("--end"):
			end=int(value)
		elif opt in ("-h"):
			usage()
			sys.exit()
	getUniqBlock(seq1,seq2,TE,name,start,end)
#	return seq1,seq2,stat,name
def getHORBlock(pwd1,pwd2,pwd3,name,start,end):
	#pwd1,pwd2,pwd3,name=getoptions()
	pwdo=name+'.block.fa'
	pwdo1=name+'.block.units.fa'
	pwdo2=name+'.block.stat.txt'
	fileo=open(pwdo,'w')
	fileo1=open(pwdo1,'w')
	fileo2=open(pwdo2,'w')
	file1=open(pwd3,'r')
	list1=[]
	for row in file1:
		row1=row.rstrip().split(' ')
		row2=[]
		for r in row1:
			if r!='':
				row2.append(r)
		list1.append([int(row2[5]),int(row2[6])])
	list21=sorted(list1,key=(lambda x:x[0]))
		
	file1.close()
	#print (len(list1))
	dict2={}
	dict3={}
	dict4={}
	for gseq in SeqIO.parse(pwd2,'fasta'):
		ID=gseq.id.split('.unit-')[0]
		loc=gseq.description.split('\t')[-2].split('-')
		if ID in dict2:
			list1t=dict2[ID]
			list2t=dict3[ID]
		else:
			list1t=[]
			list2t=[]
		list1t.append(int(loc[0])+start-1)
		list1t.append(int(loc[1])+start-1)
		list2t.append([gseq.id,int(loc[0])+start-1,int(loc[1])+start-1])
		dict2[ID]=list1t
		dict3[ID]=list2t
		dict4[gseq.id]=gseq
	list11=[]
	print (len(dict2))
	block={}
	for gseq in SeqIO.parse(pwd1,'fasta'):
		print (gseq.id)
		if gseq.id not in dict2:
			continue
		i=gseq.description.split('\t')
		s=min(dict2[gseq.id])
		e=max(dict2[gseq.id])
		cp=dict3[gseq.id]
		#print (s,e,cp)
		list1t=[]
		for i in list21:
			s1,e1=i
			x=range(s,e+1,1)
			y=range(s1,e1+1,1)
			t=set(x).intersection(set(y))
			if len(t)!=0:
				list1t.append(i)
		list2t=sorted(list1t,key=(lambda x:x[0]))
		#print (len(list2t))
		list3t=[]	
		s0=s
		for i in list2t:
			s1,e1=i
			if s1>s0:
				list3t.append([s0,s1])
			s0=e1
		if e>s0:
			list3t.append([s0,e])
	#	print(s,e,gseq.id,len(cp))
		for i in list3t:
			s1,e1=i
			list4t=[]
			for i1 in cp:
				n,s0,e0=i1
				x=range(s1,e1+1,1)
				y=range(s0,e0+1,1)
				t=set(x).intersection(set(y))
				if len(t)!=0:
					list4t.append(i1)
	#		print (len(list4t))
			list11.append([gseq.id,s1,e1,list4t,gseq])
	#	print(s,e,list3t,gseq.id,len(cp))
	list2=sorted(list11,key=(lambda x:x[2]))
	#print (len(list2))
	l1=list2[0]
	ID,s,e,cp,Seq=l1
	l1=[[ID],s,e,[cp],[Seq]]
	lenall=0
	l0=[]
	for l2 in list2[1:]:
		IDlist,s,e,cplist,Seqlist=l1
		ID1,s1,e1,cp1,Seq1=l2
		x=range(s,e+1,1)
		y=range(s1,e1+1,1)
		t=set(x).intersection(set(y))
		if len(t)!=0:
			IDlist.append(ID1)
			cplist.append(cp1)
			Seqlist.append(Seq1)
			
			l1=[IDlist,min([s,s1,e,e1]),max([s,s1,e,e1]),cplist,Seqlist]
		#	if cp1>cp:
		#		l1=l2
		#	elif cp1==cp and len(Seq1.seq)>len(Seq.seq):
		#		l1=l2
		else:
			l0.append(l1)
			l1=[[ID1],s1,e1,[cp1],[Seq1]]
	l0.append(l1)
	print ("Total Block:",len(l0))
	b0=1
	for i in l0:
		IDlist,s,e,cplist,Seqlist=i
		#if cp==[]:
		#	continue
		cpall=[]
		for cp in cplist:
			cpall=cpall+cp
			print (cp)
		print (cplist)
		if cpall==[]:
			continue
		for Seq in Seqlist:
			Seqo=Seq.description.split('\t')
			Seqo[2]=str(s)
			Seqo[3]=str(e)
			cp1=str(len(Seq.seq))+"-"+str(len(cpall))
			Seqo[1]=cp1
			Seqo1="\t".join(Seqo)+"\tBlock-"+str(b0)
			l='>'+Seqo1+"\n"+Seq.seq+"\n"
			fileo.writelines(l)
		block[b0]=cpall
		len1=[]
		loc=[]
		n=[]
		for i1 in cpall:
			gseq=dict4[i1[0]]
			loc.append(i1[1])
			loc.append(i1[2])
			l='>'+gseq.description+"\tBlock-"+str(b0)+"\n"+gseq.seq+"\n"
			fileo1.writelines(l)
			n.append(gseq.id)
			len1.append(str(len(gseq.seq)))
		print (loc)
		l="Block-"+str(b0)+"\t"+str(len(Seq.seq))+"\t"+str(len(n))+"\t"+str(min(loc))+"\t"+str(max(loc))+"\t"+str(max(loc)-min(loc)+1)+"\t"+";".join(len1)+"\t"+";".join(n)+"\n"
		#print (l)
		fileo2.writelines(l)
		b0+=1
	fileo.close()
	fileo1.close()
	fileo2.close()
#getoptions()
