import re
import Bio
from Bio import SeqIO

def uniqStandards(pwd1,pwd2):
	file1=open(pwd2,'w')
	list1=[]
	for gseq in SeqIO.parse(pwd1,'fasta'):
		loci=gseq.description.split('\t')[1].split('-')
		loc=loci[0]+"-"+loci[1]
		if loc not in list1:
			list1.append(loc)
			l='>'+gseq.description+"\n"+gseq.seq+"\n"
			file1.writelines(l)
	file1.close()
