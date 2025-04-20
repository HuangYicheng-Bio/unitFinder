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
def getHORBlock(con,pwd1,pwd4,pwd3,name,start,end):
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
		if "(0)" not in row:
			continue
		if 'Simple_repeat' in row or 'Low_complexity' in row:
			continue
		row2=[]
		for r in row1:
			if r!='':
				row2.append(r)
		list1.append([int(row2[5]),int(row2[6])])
		#list1=list1+list(range(int(row2[5]),int(row2[6])+1,1))
	list21=sorted(list1,key=(lambda x:x[0]))
	s0=start
	list22=[]
	for i in list21:
		s1,e1=i
		if s1>s0:
			list22.append([s0,s1])
		s0=e1
	if s0<end:
		list22.append([s0,end])
	file1.close()
	print (len(list22))
	dict1={}

	for gseq in SeqIO.parse(pwd1,'fasta'):
		loc=gseq.description.split('\t')[-2].split('-')
		ID=gseq.id.split('.unit-')[0]
		s=int(loc[0])+start-1
		e=int(loc[1])+start-1

		for i in list22:
			s1,e1=i
			#x=range(s,e+1,1)
			#y=range(s1,e1+1,1)
			#t=set(x).intersection(set(y))
			t=0
			if s1>=s and s1<=e:
				t=1
			if e1>=s and e1<=e:
				t=1
			if s>=s1 and s<=e1:
				t=1
			if e>=s1 and e<=e1:
				t=1
			if t!=0:
				k=str(s1)+"-"+str(e1)
				if k in dict1:
					list1t=dict1[k]
				else:
					list1t=[[],[],[],[]]
				list1t[0].append(s)
				list1t[0].append(e)
				list1t[1].append(gseq.id)
				list1t[2].append(gseq)
				if ID not in list1t[3]:
					list1t[3].append(ID)
				dict1[k]=list1t
	print ("CentGm-1",len(dict1))
	for k,v in dict1.items():
		print (len(v[1]))
	dict4={}
	for gseq in SeqIO.parse(pwd4,'fasta'):
		loc=gseq.description.split('\t')[-2].split('-')
		ID=gseq.id.split('.unit-')[0]
		s=int(loc[0])+start-1
		e=int(loc[1])+start-1
		
		for i in list22:
			s1,e1=i
			t=0
			#x=range(s,e+1,1)
			#y=range(s1,e1+1,1)
			#t=set(x).intersection(set(y))
			if s1>=s and s1<=e:
				t=1
			if e1>=s and e1<=e:
				t=1
			if s>=s1 and s<=e1:
				t=1
			if e>=s1 and e<=e1:
				t=1
			if t!=0:
				k=str(s1)+"-"+str(e1)
				if k in dict4:
					list1t=dict4[k]
				else:
					list1t=[[],[],[],[]]
				list1t[0].append(s)
				list1t[0].append(e)
				list1t[1].append(gseq.id)
				list1t[2].append(gseq)
				if ID not in list1t[3]:
					list1t[3].append(ID)
				dict4[k]=list1t
	print ("CentGm-4",len(dict4))
	for k,v in dict4.items():
		print (len(v[1]))
	l0=[]
	
	for i in list22:
		s1,e1=i
		k=str(s1)+"-"+str(e1)
		if k in dict1:
			b1=dict1[k]
		else:
			b1=[]
		if k in dict4:
			b4=dict4[k]
		else:
			b4=[]
		if b1!=[] and b4!=[]:
			if min(b1[0])>min(b4[0]):
				l0.append([b4[1],min(b4[0]),max(b4[0]),len(b4[1]),b4[2],b4[3],'CentGm-4'])
				l0.append([b1[1],min(b1[0]),max(b1[0]),len(b1[1]),b1[2],b1[3],'CentGm-1'])
			else:
				l0.append([b1[1],min(b1[0]),max(b1[0]),len(b1[1]),b1[2],b1[3],'CentGm-1'])
				l0.append([b4[1],min(b4[0]),max(b4[0]),len(b4[1]),b4[2],b4[3],'CentGm-4'])
		else:
			if b1!=[]:
				l0.append([b1[1],min(b1[0]),max(b1[0]),len(b1[1]),b1[2],b1[3],'CentGm-1'])
			if b4!=[]:
				l0.append([b4[1],min(b4[0]),max(b4[0]),len(b4[1]),b4[2],b4[3],'CentGm-4'])
	print ("Total Block:",len(l0))
	b0=1
	for i in l0:
		IDlist,s,e,cpnumber,Seqlist,Seqlistcon,type1=i
		
		for Seq in SeqIO.parse(con,'fasta'):
			if Seq.id in Seqlistcon:
				Seqo=Seq.description.split('\t')
				Seqo[2]=str(s)
				Seqo[3]=str(e)
				cp1=str(len(Seq.seq))+"-"+str(cpnumber)
				Seqo[1]=cp1
				Seqo1="\t".join(Seqo)+"\tBlock-"+str(b0)+"\t"+type1
				l='>'+Seqo1+"\n"+Seq.seq+"\n"
				fileo.writelines(l)
		
		len1=[]
		loc=[]
		n=[]
		for gseq in Seqlist:
			loc1=gseq.description.split('\t')[-2].split('-')
			loc.append(int(loc1[0]))
			loc.append(int(loc1[1]))
			l='>'+gseq.id+"\t"+type1+"\t"+"\t".join(gseq.description.split("\t")[1:])+"\tBlock-"+str(b0)+"\n"+gseq.seq+"\n"
			#l='>'+gseq.description+"\tBlock-"+str(b0)+"\t"+type1+"\n"+gseq.seq+"\n"
			fileo1.writelines(l)
			n.append(gseq.id)
			len1.append(str(len(gseq.seq)))
		print (min(loc),max(loc))
		l="Block-"+str(b0)+"\t"+type1+"\t"+str(len(Seq.seq))+"\t"+str(len(n))+"\t"+str(min(loc))+"\t"+str(max(loc))+"\t"+str(max(loc)-min(loc)+1)+"\t"+";".join(len1)+"\t"+";".join(n)+"\n"
		#print (l)
		fileo2.writelines(l)
		b0+=1
	fileo.close()
	fileo1.close()
	fileo2.close()
#getoptions()
