import re
import sys
import Bio
from Bio import SeqIO
from matplotlib import colors
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import getopt
import sys
def usage():
	print ("--in *.loci.fa")
	print ("--seq *.fa")
	print ("--txt *.block.txt")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","seq=","txt="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	seq3=''
	name=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		if opt in ("--txt"):
			seq3=value
		if opt in ("--seq"):
			seq2=value
		elif opt in ("--name"):
			name=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,seq2,seq3,name

pwd1,pwd2,pwd3,name=getoptions()

litt=[name]
#litt=['Chr08']#,'Chr08','Chr03','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19','Chr20',]
#pwd1='Chr08-29459939-32372411.uniq.fa'
#pwd1=' Chr01.raw.loci.fa'
pwdo=name+".selected.element.txt"
#pwd2='../../../1-Loci/Chr08/Chr08-29459939-32372411.fa'
#pwd3='Chr08-29459939-32372411.uniq.txt'
fileo=open(pwdo,'w')

def heatmap(data,ax=None,cbar_kw={}, cbarlabel="",**kwargs):
	if not ax:
		ax = plt.gca()
	im = ax.imshow(data, **kwargs)
	cbar = ax.figure.colorbar(im,shrink=0.5, norm=colors.LogNorm(vmin=0, vmax=500))#,fraction=0.05)
	cbar.ax.set_ylabel(cbarlabel,fontsize =60,rotation=-90, va="bottom")
	cbar.ax.yaxis.set_tick_params(labelsize=40)
	ax.tick_params(top=False, bottom=True,labeltop=False, labelbottom=True)
	plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")
	#ax.spines[:].set_visible(False)
	#ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
	#ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
	ax.xaxis.set_tick_params(labelsize=40)
	ax.yaxis.set_tick_params(labelsize=40)
	ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
	#ax.tick_params(which="minor", top=False,bottom=True, left=False)
	return im, cbar


nn=0
mn=0
for Chrn in litt:
	fig,ax=plt.subplots(figsize=(82.7,117))
	for gseq in SeqIO.parse(pwd2,'fasta'):
		ChrnLen=len(gseq.seq)
		Chrn1=gseq.id
	n=0
	heat=[]
	score=[]
	
	xlable=[]
	while n<=ChrnLen:
		n1=[0.0]*1000
		n2=[-100]*1000
		heat.append(n1)
		score.append(n2)

		xlable.append(str(n+1)+"-"+str(n+1001))
		n+=1000

	print (len(heat),len(heat)*len(n1),Chrn,ChrnLen)
	n=0
	for gseq1 in SeqIO.parse(pwd1,'fasta'):
		row1=gseq1.description.split('\t')
		s=int(row1[2])
		e=int(row1[3])
		f=0
		filet=open(pwd3,'r')
		flag=0
		ID1=gseq1.id.split('_Repeat')
		ID2='_Repeat'.join(ID1[:2])
		for ru in filet:
			ru1=ru.rstrip().split('\t')
			IDlist=ru1[-1].split(';')
			if len(IDlist)!=1:
				for iID in IDlist:
					if ID2==iID:
						flag=1
						print (ru)
						len1=int(ru1[1])
		filet.close()
		if flag==0:
			continue
		for i in range(s,e+1):
			y=int(i/1000)
			x=int(i-y*1000)
			
			if score[y-1][x-1]==-100.0:
				heat[y-1][x-1]=len1
				score[y-1][x-1]=float(row1[9])
				f=1
				
			else:
				sc=score[y-1][x-1]
				if sc>float(row1[2]):
					heat[y-1][x-1]=len1
					score[y-1][x-1]=float(row1[9])
					f=1
		if f==1:
			l='>'+gseq1.description+"\t"+str(len(gseq.seq))+"\n"
			#l='>'+gseq1.description+"\n"+gseq.seq+"\n"
			print (gseq1.description)
			fileo.writelines(l)
			print (gseq1.description)

	
#	norm=mpl.colors.BoundaryNorm(np.linspace(-500, 500), 1000)
	im,cbar = heatmap(np.array(heat),ax=ax,cmap="ocean_r",cbarlabel="Length of Tandem Repeats")
	#im,cbar = heatmap(heat1,ax=ax,cmap="tab20c_r",cbarlabel="Length of CentGm-4")
	#ax.set_xlabel(Chrn1,fontsize=60)
	ax.set_title('>'+Chrn1,fontsize=100,loc='left')
	
	figname=Chrn+"-CentGm.pdf"
	fig.savefig(figname)
	plt.close()
print ("welldone")
fileo.close()

