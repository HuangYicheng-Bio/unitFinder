import re
import os
import Bio
from Bio import SeqIO
from mummer_type import mummer
import getopt
import sys
def usage():
	print ("--in *.fa")
	print ("--out *.group.txt")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","out=","name="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	name=''
	trf=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		elif opt in ("--out"):
			seq2=value
		elif opt in ("--name"):
			name=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	typeMerge(seq1,seq2,name)
	return seq1,seq2,name

def typeMerge(pwd1,pwd2,name):
#pwd1,pwd2,name=getoptions()

	file2=open(pwd2,'w')
	dist1={}
	list1=[]
	for gseq in SeqIO.parse(pwd1,'fasta'):
		if len(gseq.seq)>10:
			list1.append(gseq.id)
			dist1[gseq.id]=gseq

	seq1=list1[0]
	dictg={}
	n=[]
	listg=[]
	for seq1 in list1:
		if seq1 in n:
			continue
		else:
			print (len(n),float(len(n))/len(list1))
		dictg[seq1]=[]
		n.append(seq1)
		listt=[seq1]
		j=0
		n1=[]
		while j<len(list1):
			seq2=list1[j]
		#for seq2 in list1:
			tag=0
			for v in listt:
				n2=v+"-"+seq2
				if n2 in n1:
					continue
				else:
					n1.append(n2)
				if v!=seq2 and seq2 not in n:
					tag=0
					gseq1=dist1[v]
					gseq2=dist1[seq2]
					filen1=name+'_type/'
					pwda="typetemp1.fasta"
					filet=open(pwda,'w')
					l='>'+gseq1.id+"\n"+gseq1.seq+"\n"
					filet.writelines(l)
					filet.close()
					pwdb="typetemp2.fasta"
					filet=open(pwdb,'w')
					l='>'+gseq2.id+"\n"+gseq2.seq+"\n"
					filet.writelines(l)
					filet.close()
			
					mummer(pwda,pwdb,'typetemp3')
					match1=[]
					match2=[]
					file1=open('typetemp3.delta.filter.coords','r')
					for row in file1:
						row1=row.rstrip().split('\t')
						if float(row1[6])>=80:#80.0:
							s=int(row1[2])-1
							e=int(row1[3])
							if s>e:
								s=int(row1[3])-1
								e=int(row1[2])
							match1+=range(s,e+1,1)
							s=int(row1[0])-1
							e=int(row1[1])
							match2+=range(s,e+1,1)
						if len(set(match1))/float(row1[8])>=0.8 or len(set(match2))/float(row1[7])>=0.8:
							tag=1
							break
					file1.close()
					commandline="rm typetemp*"
					os.system(commandline)
					#if v=='CP126426.1_Chr01-14942553-16318859_Repeat_1-130' or seq2=='CP126442.1_Chr17-26178551-28319679_Repeat_1-315':
					#	print (v,seq2,tag,len(set(match1)),len(set(match2)))
					#if v=='CP126426.1_Chr01-14942553-16318859_Repeat_1-130' and  seq2=='CP126442.1_Chr17-26178551-28319679_Repeat_1-315':
					
					#	sys.exit()
					if tag==1:
						listt.append(seq2)
						n.append(seq2)
						j=0
						print (v,seq2,tag,len(set(match1)),len(set(match2)),float(row1[7]),float(row1[8]))
						break
			if tag==0:
				j+=1
		listg.append(listt)
	j=1
	os.system('mkdir '+name+"-TypeGroup")
	#print ('Type number:',len(dictg))
	pwdfile=name+"-TypeGroup/"
	for v in listg:
		filet=open(pwdfile+"Type-"+str(j)+"-all.fasta",'w')
		lc=[]
		for i in list(set(v)):
			gseq=dist1[i]
			l1='>'+gseq.description+"\tType-"+str(j)+'\n'+gseq.seq+'\n'
			filet.writelines(l1)
			cp=float(gseq.description.split('\t')[2].split('-')[1])
			if lc==[]:
				lc=[cp,gseq]
			else:
				if cp>lc[0]:
					lc=[cp,gseq]
		filet.close()
		
		filet=open(pwdfile+"Type-"+str(j)+"-typic.fasta",'w')
		gseqt=lc[1]
		l='>'+gseqt.description+"\tType-"+str(j)+'\n'+gseqt.seq+'\n'
		filet.writelines(l)
		filet.close()
		l='Type-'+str(j)+'\t'+gseqt.id+"\t"+str(len(gseqt.seq))+"\t"+str(int(float(lc[0])))+"\t"+';'.join(v)+"\n"
		file2.writelines(l)
		j+=1
	file2.close()
getoptions()
