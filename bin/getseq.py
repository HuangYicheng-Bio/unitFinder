import re
import os
import sys
import getopt
import Bio
from Bio import SeqIO

def usage():
	print ("--input| -i\tSequence file")
	print ("--ID | -n\t Seqs_ID_file | ID")
	print ("--out | -o\toutfile")
	print ("--start | -s start")
	print ("--end | -e end")
	print ("--flag | -t\tstrand")

def getoptions():
	try:
		opts,args=getopt.getopt(sys.argv[1:],"I:n:o:s:e:o:t:h",["input=","ID=","out=","start=","end=","flag=","help"])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()
	input1=""
	ID=""
	out=""
	start=-1
	end=-1
	flag="+"
	
	if len(opts)==0:
		print ("NO parameter!")
		print ("Use '--help | -h for same information")
		sys.exit()
	
	for opt,value in opts:
		if opt in ("--input","-i"):
			input1=value
		elif opt in ("--ID","-n"):
			ID=value
		elif opt in ("--start","-s"):
			start=int(value)
		elif opt in ("--end","-e"):
			end=int(value)
		elif opt in ("--flag","-t"):
			flag=value
		elif opt in ("--out","-o"):
			out=value
		elif opt in ("--help","-h"):
			usage()
			sys.exit()
	return input1,ID,out,start-1,end,flag


def getseq(input1,start,end,out):
	start=start-1
	fileo=open(out,'w')
	for gseq in SeqIO.parse(input1,'fasta'):
		Seqo=gseq.seq[start:end]
		l=">"+gseq.id+"-"+str(start)+"-"+str(end)+"\n"+str(Seqo)+"\n"
		fileo.writelines(l)
	fileo.close()
	print ("welldone")
						
