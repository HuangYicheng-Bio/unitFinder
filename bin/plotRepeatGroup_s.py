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
	print ("--group *.grouped.txt")
	print ("--name ID")
	print ("--start startloci")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","start=","group="])
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
		elif opt in ("--group"):
			seq2=value
		elif opt in ("--start"):
			start=int(value)
		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,seq2,start,name

pwd1,pwd2,start,name=getoptions()
fig,ax=plt.subplots(figsize=(10,5))

dict1={}
file1=open(pwd2,'r')
ty1=[]
for row in file1:
	ID=row.rstrip().split('\t')
	list1=ID[-1].split(';')
	dict1[ID[0]]=list1
	ty1=ty1+list1	
file1.close()

norm=colors.Normalize(1, len(dict1))
cmap = plt.colormaps["viridis"]

#print (min(lc1), max(lc1))
dictc={}
for gseq in SeqIO.parse(pwd1,'fasta'):
	loci=gseq.description.split('\t')[-1].split('-')
	s=start+int(loci[0])-1
	e=start+int(loci[1])-1
	t=int(len(gseq.seq))
	x=[s,e]
	y=[t,t]
	ID=gseq.id.split('.unit-')[0]
	type1=''
	for k,v in dict1.items():
		if ID in v:
			type1=k
			break
	if type1!='':
		t=int(type1.split('-')[1])
		c= plt.cm.viridis(norm(t))
		ax.plot(x,y,color='black')
		dictc[gseq.id]=c
ax.set_xlabel("locition")
ax.set_ylabel("Units Length")

for gseq in SeqIO.parse(pwd1,'fasta'):
	if gseq.id in dictc:
		loci=gseq.description.split('\t')[-1].split('-')
		s=start+int(loci[0])-1
		e=start+int(loci[1])-1
		c=len(gseq.seq)
		vertices = [(s, c-20), (s, c+20), (e, c+20), (e, c-20), (0, 0)]
		codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
		path = Path(vertices, codes)
		c=dictc[gseq.id]
		pathpatch = PathPatch(path, alpha=0.5,facecolor='none', edgecolor=c)
	
		#ax.add_patch(pathpatch)
		ax.axvspan(s, e, facecolor=c,alpha=0.1)

fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),label='Type',ax=ax)
fig.suptitle(name+" (#Type: "+str(len(dict1))+")", size=14)
plt.tight_layout()
figname=name+".pdf"
fig.savefig(figname)
plt.close()
print ("welldone")
