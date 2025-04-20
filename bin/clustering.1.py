import re
import numpy as np
from numpy import unique
from numpy import where
from sklearn.cluster import AffinityPropagation
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import statistics
import getopt
import sys
def usage():
	print ("--in *.delta.coords")
	print ("--RP *.fa")
	print ("--fa *.fa")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","fa","--RP=","name="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	name=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		elif opt in ("--fa"):
			seq2=value

		elif opt in ("--name"):
			name=value

		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,seq2,name

pwd1,pwd2,name=getoptions()

file1=open(pwd1,'r')
file2=open(pwd2,'r')
list1=[]
list2=[]
i0=0
x0=[]
y0=[]
x1=[]
y1=[]


for row in file1:
	row1=row.rstrip().split('\t')
	s=int(row1[0])
	e=int(row1[1])
	s1=int(row1[2])
	e1=int(row1[3])
	m1=int((s+e)/2.0)
	m2=int((s1+e1)/2.0)
	x1.append(m1)
	y1.append(m2)
l1=[]
l1s=[]
l1e=[]
fig, ax = plt.subplots(1, 3,sharey=True,sharex=True,figsize=(45, 15))
j0=0
s=[]
cl=[]
x=[]
y=[]
sl=[]
cl=[]
for cluster in clusters:
	row_ix = where(yhat == cluster)
	l2=[]
	l2s=[]
	l2e=[]
	xt=[]
	yt=[]
	for i in row_ix[0]:
		l2.append(np.mean(list1[i]))
		l2s.append(list2[i][0])
		l2e.append(list2[i][1])
		x.append(list1[i][0])
		xt.append(list1[i][0])
		y.append(list1[i][1])
		yt.append(list
scatter = ax[2].scatter(x, y,c=sl,edgecolors='none')

axins1 = inset_axes(ax[2],width="3%",height="50%",loc="upper left")
ax[1].legend(loc="upper left")

cbar = fig.colorbar(scatter,cax=axins1)
cbar.ax.set_ylabel('cluster number')
ax[0].scatter(x1,y1,s=0.5)

i0=0
l=[]
for i in l1:
	lo=str(np.mean(i))+"\t"+str(min(i))+"\t"+str(max(i))+"\t"+str(max(i)-min(i))+"\t"+str(len(i))+"\t"+str(l1s[i0])+"\t"+str(l1e[i0])+"\n"
	if l==[]:
		l=[len(i),l1s[i0],l1e[i0]]
	else:
		if len(i)>l[0]:
			l=[len(i),l1s[i0],l1e[i0]]
	file2.writelines(lo)
	i0+=1

codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
x=l[1]
y=l[2]
vertices = [(x, x), (x, y), (y, y), (y, x), (0, 0)]

path = Path(vertices, codes)
pathpatch = PathPatch(path, facecolor='none', edgecolor='Red')
ax[0].add_patch(pathpatch)
ax[0].grid(True)
ax[0].set_title(name+' Original Alignment')

path = Path(vertices, codes)
pathpatch = PathPatch(path, facecolor='none', edgecolor='Red')

ax[0].xaxis.set_major_formatter(lambda x, pos: str(x/1000000)+" Mbp")
ax[0].yaxis.set_major_formatter(lambda x, pos: str(x/1000000)+" Mbp")

figname=name+".png"
plt.subplots_adjust(hspace=0.01,wspace=0.01)
fig.savefig(figname, bbox_inches='tight')
plt.close()
print ('Potential Centromeric region ',str(l[1]),str(l[2]))
file1.close()
file2.close()
print ('welldone')
