import re
import numpy as np
import Bio
from Bio import SeqIO
import re
import getopt
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
pwd1='CentGM4.fa'
pwd2='CentGm-4.type.fa'
file1=open(pwd2,'w')
list1=[]
j=0
listc=[]
for gseq in SeqIO.parse(pwd1,'fasta'):
	Seq1=list(str(gseq.seq))
	if list1==[]:
		list1=Seq1
	#	print (list1)
	#	sys.exit()
	else:
		list2=zip(list1,Seq1)
		list1=[]
		for i in list2:
			i1=list(eval(str(i).replace("[", '').replace("]", '')))
			list1.append(i1)
			if j==1:
				listc.append([])
		#print (list1)
	j+=1

print (listc)
for gseq in SeqIO.parse(pwd1,'fasta'):
	Seq1=list(str(gseq.seq))
	j=0
	while j<len(Seq1):
		print (listc[j],Seq1[j])
		if Seq1[j]!='-':
			listc[j].append(gseq.id.split('-')[0])
		j+=1
print (listc)

listq=[]
j=0
category_names = ['A','T','G','C']
results={}
for i in list1:
	lt=[]
	for i1 in i:
		if i1!='-':
			lt.append(i1.upper())
	m=[]
	m1=[]
	jl=[]
	for i2 in ['A','T','G','C']:
		i3=lt.count(i2)
		jl.append(i3)
		if m==[]:
			m=[i2,i3]
		else:
			if m[1]<i3:
				m=[i2,i3]
			#elif m[1]==i3:
			#	print (m,i2,i3)
			#	sys.exit()
		m1.append([i2,i3])		
	listq.append(m[0])
	results[j]=jl
	print (m1,m[0],m[1],len(set(listc[j])),set(listc[j]))
	j+=1
SeqO=''.join(listq)
print (len(SeqO))
for gseq in SeqIO.parse(pwd1,'fasta'):
	print (len(gseq.seq))
l='>CentGm-4\n'+SeqO+"\n"
file1.writelines(l)
file1.close()

fig, ax = plt.subplots(1, 1, figsize=(12, 14))
labels = list(results.keys())
data = np.array(list(results.values()))
data_cum = data.cumsum(axis=1)
#category_colors = plt.colormaps['RdYlGn'](np.linspace(0.15, 0.85, data.shape[1]))
category_colors =["red","green","orange","blue"]
fig, ax = plt.subplots(figsize=(5,10))
ax.invert_yaxis()
ax.xaxis.set_visible(False)
ax.set_xlim(0, np.sum(data, axis=1).max())
for i, (colname, color) in enumerate(zip(category_names, category_colors)):
	widths = data[:, i]
	starts = data_cum[:, i] - widths
	rects = ax.barh(labels, widths, left=starts, height=0.5,label=colname, color=color)
#	r, g, b, _ = color
#	text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
#	ax.bar_label(re
ax.legend(ncols=len(category_names), bbox_to_anchor=(0, 1),loc='lower left', fontsize='small')
figname="CentGm-4.unify.pdf"
fig.savefig(figname)
print ("welldone")

