import re
import getopt
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import Bio
from Bio import SeqIO
def usage():
	print ("--in *.units.fa")
	print ("--all *Consensus.fa")
	print ("--name ID")
	print ("--start startloci")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","start=","all="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
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
		elif opt in ("--all"):
			seq2=value
		elif opt in ("--start"):
			start=int(value)
		elif opt in ("-h"):
			usage()
			sys.exit()
	plotRepeatBlock(seq1,seq2,name,start)
	#return seq1,seq2,start,name

def plotRepeatBlock(pwd1,pwd2,name,start):
	#pwd1,pwd2,start,name=getoptions()
	file0=open(name+'.stat.txt','w')
	fig,ax=plt.subplots(figsize=(10,5))

	dict1={}

	for gseq in SeqIO.parse(pwd1,'fasta'):
		ID=gseq.id.split('.unit-')[0]
		if ID in dict1:
			list1=dict1[ID]
		else:
			list1=[]
		#if gseq.id not in list1:
		list1.append(gseq.id)
		dict1[ID]=list1
		
	
	lc1=[]
	for k,v in dict1.items():
		if len(v)!=len(set(v)):
			print (v)
			print ('error')
			sys.exit()
		l=k+'\t'+str(len(v))+"\t"+";".join(v)+"\n"
		file0.writelines(l)
		if v not in lc1:
			lc1.append(len(v))
	
	file0.close()
	
	norm=colors.Normalize(min(lc1), max(lc1))
	cmap = plt.colormaps["viridis"]
	
	print (min(lc1), max(lc1))
	for gseq in SeqIO.parse(pwd1,'fasta'):
		loci=gseq.description.split('\t')[-1].split('-')
		s=start+int(loci[0])-1
		e=start+int(loci[1])-1
		t=int(len(gseq.seq))
		x=[s,e]
		y=[t,t]
		ID=gseq.id.split('.unit-')[0]
		n=dict1[ID]
		c= plt.cm.viridis(norm(len(n)))
		ax.plot(x,y,color=c)
	ax.set_xlabel("locition")
	ax.set_ylabel("Units Length")

	for gseq in SeqIO.parse(pwd2,'fasta'):
		loci=gseq.description.split('\t')
		s=start+int(loci[2])-1
		e=start+int(loci[3])-1
		c=len(gseq.seq)
		vertices = [(s, c-20), (s, c+20), (e, c+20), (e, c-20), (0, 0)]
		codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
		path = Path(vertices, codes)
		
		pathpatch = PathPatch(path, facecolor='none', edgecolor='Red')
	
		ax.add_patch(pathpatch)
		#ax.axvspan(s, e, facecolor='#2ca02c')

	fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),label='Repeat Number',ax=ax)
	fig.suptitle(name+" (#Repeat: "+str(len(dict1))+")", size=14)
	plt.tight_layout()
	figname=name+".pdf"
	fig.savefig(figname)
	plt.close()
	print ("welldone")
#getoptions()
