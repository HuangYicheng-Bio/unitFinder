import re
import Bio
from Bio import SeqIO
import sys
import getopt
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def usage():
	print ("--in\t*.fasta")
	print ("--out|-o\toutput")
	print ("-h\tDisplay help")


def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "o:h",["in=","out="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	if len(opts)==0:
		print ("PARAMETER ERROR!")
		print ("Use '--h' for same information")
		sys.exit()

	fastafile=""
	out=""

	for opt,value in opts:
		if opt in ("--in"):
			fastafile=value
		elif opt in ("--out","-o"):
			out=value

		elif opt in ("-h"):
			usage()
			sys.exit()

	return fastafile,out

def length(fastapwd,outf):

	i=0
	l=[]
	for gseq in SeqIO.parse(fastapwd,'fasta'):
		i=i+len(gseq.seq)
		l.append(len(gseq.seq))
	print ("Total Length",i)
	print ("welldone")


	fig, ax = plt.subplots()
	figname1=outf+".pdf"
	ax.hist(l,bins=500)
	plt.title(outf)
	fig.savefig(figname1)
	plt.close()
