import re
import os
import Bio
from Bio import SeqIO
from mummer_type import mummer
import getopt
import sys
def usage():
	print ("--in *.fa")
	print ("--out *.coords")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","out="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		elif opt in ("--out"):
			seq2=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,seq2

pwd1,pwd2=getoptions()

file2=open(pwd2,'w')
dist1={}
list1=[]
for gseq in SeqIO.parse(pwd1,'fasta'):
	list1.append(gseq.id)
	dist1[gseq.id]=gseq

n=[]

for seq1 in list1:
	n.append(seq1)
	for seq2 in list1:
		if seq2 not in n:				
			gseq1=dist1[seq1]
			gseq2=dist1[seq2]
			pwda="temp1.fasta"
			filet=open(pwda,'w')
			l='>'+gseq1.id+"\n"+gseq1.seq+"\n"
			filet.writelines(l)
			filet.close()
			pwdb="temp2.fasta"
			filet=open(pwdb,'w')
			l='>'+gseq2.id+"\n"+gseq2.seq+"\n"
			filet.writelines(l)
			filet.close()
	
			mummer(pwda,pwdb,'temp3')
			file1=open('temp3.delta.filter.coords','r')
			for row in file1:
				file2.writelines(row)

			file1.close()
			commandline="rm temp*"
			os.system(commandline)	
file2.close()
