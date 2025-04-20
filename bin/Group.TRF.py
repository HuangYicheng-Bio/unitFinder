import re
import numpy as np
from numpy import unique
from numpy import where
from sklearn.cluster import AffinityPropagation,DBSCAN
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
	print ("--trf *.dat")
	print ("--out *.txt")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","out=","name=","trf="])
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
		elif opt in ("--trf"):
			trf=value
		elif opt in ("--name"):
			name=value

		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,seq2,name,trf

pwd1,pwd2,name,trf=getoptions()

file1=open(trf,'r')
pwdtrf=name+'.TRF.txt'
file1o=open(name+'.TRF.txt','w')

for row in file1:
	row1=row.rstrip().split(' ')
	if len(row1)>=12:
		a=row1[0]
		b=row1[1]
		cpnumber=int(float(row1[3]))
		l=a+"\t"+b+"\t"+a+"\t"+b+"\t"+row1[2]+"\t"+str(cpnumber)+"\n"
		file1o.writelines(l)
file1.close()
file1o.close()

file1=open(pwdtrf,'r')
cn=[]
dict1={}
for row in file1:
	row1=row.rstrip().split('\t')
	if int(row1[-1])<100:
		continue
	cn.append(int(row1[-1]))
	if int(row1[-1]) in dict1:
		lt=dict1[int(row1[-1])]
	else:
		lt=[]
	lt.append([int(row1[0]),int(row1[1])])
	dict1[int(row1[-1])]=lt
file1.close()

xA=np.array(cn)
mean1=np.mean(xA)
std1=np.std(xA)
out1=[]
m=np.percentile(xA,(25,75))
IQR=m[1]-m[0]
t=3*IQR
t0=min(max(xA),m[1]+t)
for i in xA:
	if i>t0:
		out1.append(i)
out1=sorted(list(set(out1)))

lt1=[]
for i in out1:
	lt=dict1[i]
#	print (dict1[i])
	for j in lt:
		j.append(i)
	lt1=lt1+lt
lTRF=[]
print (lt1)
l1=sorted(lt1,key=(lambda x:x[0]))
for l1t in l1:
	a1,b1,c1=l1t
	if lTRF==[]:
		lTRF=[[a1,b1,c1]]
	else:
		ltn=0
		flag=0
		while ltn<len(lTRF):
			a,b,c=lTRF[ltn]
			xt1=range(a1-1000000,b1+1000000,1)
			yt1=range(a-1000000,b+1000000,1)
			zt1=list(set(xt1).intersection(set(yt1)))
			if len(zt1)!=0:
				a2=min(a1,a)
				b2=max(b1,b)
				c2=c1+c
				lTRF[ltn]=[a2,b2,c2]
				flag=1
			#	print ([a2,b2,c2])
				break
			ltn+=1
		if flag==0:
			lTRF.append([a1,b1,c1])
for i in lTRF:
	a,b,c=i
	#print (a,b,b-a,c)

file1=open(pwd1,'r')
file2=open(pwd2,'w')
list1=[]
list2=[]

x0=[]
y0=[]

for row in file1:
	row1=row.rstrip().split('\t')
	
	s=int(row1[0])
	e=int(row1[1])
	s1=int(row1[2])
	e1=int(row1[3])
	m1=int((s+e)/2.0)
	m2=int((s1+e1)/2.0)
	
	n01=min([s,e,s1,e1])
	n02=max([s,e,s1,e1])
	list1.append([m1,m2])
	x0.append(m1)
	y0.append(m2)
	list2.append([n01,n02])
	
#print (len(list1))

l1=[]
l1s=[]
l1e=[]
print ('cluster number:',len(lTRF))
fig, ax = plt.subplots(1, 1,figsize=(10, 10))

ax.scatter(x0,y0,s=0.5)

l=[]
j0=1
lmax=[]
for i in lTRF:
	a,b,c=i
	lo='Cluster-'+str(j0)+"\t"+str(a)+"\t"+str(b)+"\t"+str(b-a)+"\t"+str(c)+"\n"
	l.append([c,a,b,[str(j0)]])
	if lmax==[]:
		lmax=[[c,a,b,[str(j0)]]]
	else:
		if c>lmax[0][0]:
			lmax=[[c,a,b,[str(j0)]]]
	j0+=1
	file2.writelines(lo)

lco=''
lco1=''

g=0
for l1t in l:
	codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
	x=l1t[1]
	y=l1t[2]
	lco1=lco1+'Cluster-'+" Cluster-".join(l1t[3])+" Number:"+str(l1t[0])+' Start: '+str(x)+' End: '+str(y)+"\n"
	lco=lco+'Location-'+str(g)+" #Tandem Repeats Units: "+str(l1t[0])+' Start: '+str(x)+' End: '+str(y)
	vertices = [(x, x), (x, y), (y, y), (y, x), (0, 0)]

	path = Path(vertices, codes)
	if l1t==lmax[0]:
		pathpatch = PathPatch(path, facecolor='none', edgecolor='Red')
		lco=lco+" (Potential Centromeric region)\n"
		ax.add_patch(pathpatch)
	else:
		pathpatch = PathPatch(path, facecolor='none', edgecolor='Yellow')
		lco=lco+"\n"
		#ax.add_patch(pathpatch)
	g+=1

print (lco)

ax.grid(True)
ax.set_title('Potential Centromeric region')
ax.set_xlabel(lco,loc='left')

ax.xaxis.set_major_formatter(lambda x, pos: str(x/1000000)+" Mbp")
ax.yaxis.set_major_formatter(lambda x, pos: str(x/1000000)+" Mbp")

figname=name+".png"

plt.subplots_adjust(hspace=0.01,wspace=0.01)
fig.savefig(figname, bbox_inches='tight')
plt.close()
print (lco)
file1.close()
file2.close()
print ('welldone')
