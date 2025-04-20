import Bio
import re
from Bio import SeqIO
import sys
from matplotlib import colors
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def heatmap(data,ax=None,cbar_kw={}, cbarlabel="",**kwargs):
	if not ax:
		ax = plt.gca()
	#norm=mpl.colors.BoundaryNorm(np.linspace(-0.5,5.5,7), 6)
	
	im = ax.pcolormesh(data,**kwargs)
	cbar=[]
	#cbar = ax.figure.colorbar(im,shrink=0.5,**cbar_kw)
	#cbar.ax.set_ylabel(cbarlabel,fontsize =80,rotation=-90, va="bottom")
	#cbar.ax.yaxis.set_tick_params(labelsize=60)
	ax.tick_params(top=False, bottom=True,labeltop=False, labelbottom=True)
	plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")
	#ax.spines[:].set_visible(False)
	#ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
	#ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
	ax.xaxis.set_tick_params(labelsize=20)
	ax.yaxis.set_tick_params(labelsize=20)
	ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
	ax.tick_params(which="minor", bottom=False, left=False)
	return im, cbar
def aln4(pwd1,pwd2,name):
	#pwd1='../CentGm-4-Chr02.complete.fa'
	#pwd2='CentGm-4-Chr02.type.fa'
	dict1={}
	for gseq in SeqIO.parse(pwd1,'fasta'):
		if len(gseq.seq)<=300:
			continue
		if gseq.seq in dict1:
			list1=dict1[gseq.seq]
		else:
			list1=[]
		list1.append(gseq.description)
		dict1[gseq.seq]=list1
	
	dict2={}
	num=0
	v1=[]
	list2=[]
	for k,v in dict1.items():
		list2.append([k,len(v)])
		if len(v)>=40:
			print (len(k),len(v))
			num+=len(v)
			dict2[k]=v
			for i in v:
				i1=i.split('\t')[-2].split('-')
				s=int(i1[0])
				e=int(i1[1])
				v1.append([(e+s+1)/2,i.split('\t')[0]])
	list3=sorted(list2,reverse=True,key=(lambda x:x[1]))
	print (list3[:5])
	v11=[]
	v12=[]
	dict21={}
	num=0
	num1=1
	dict4={}
	file2=open(pwd2,'w')
	for seq1 in list3[:10]:
		v=dict1[seq1[0]]
		dict21[seq1[0]]=v
		l='>'+name+'_Type-'+str(num1)+"\t"+str(len(v))+"\t"+";".join(v)+"\n"+seq1[0]+"\n"
		file2.writelines(l)
		num+=len(v)
		for i in v:
			i1=i.split('\t')[-2].split('-')
			s=int(i1[0])
			e=int(i1[1])
			v11.append([(e+s+1)/2,i.split('\t')[0]])
			v12.append([s,e,i.split('\t')[0]])
			dict4[i.split('\t')[0]]=num1
		num1+=1
	file2.close()
	#v2=sorted(v1,key=(lambda x:x[0]))
	v2=sorted(v11,key=(lambda x:x[0]))
	v21=sorted(v12,key=(lambda x:x[0]))
	print (v2[:10])
	#n0=0
	#dict3={}
	#for v3 in v2:
	#	dict3[v3[1]]=n0
	#	n0+=1
	
	loc1=[]
	v0=[v21[0][0],v21[0][1],[v21[0][-1]]]
	vabord=[]
	for v3 in v21[1:]:
		s1=v0[0]
		e1=v0[1]
		s2=v3[0]
		e2=v3[1]
		if s2-10000<=e1+10000:
			v0[1]=e2
			v0[-1].append(v3[-1])
		else:
			#loc1.append(v0)
			print (v0[0],v0[1],len(v0[-1]))
			if len(v0[-1])>10:
				#print (v0[0],v0[1],len(v0[-1]))
				loc1.append(v0)
			else:
				for va in v0[-1]:
					vabord.append(va)
			v0=[v3[0],v3[1],[v3[-1]]]
	print (v0[0],v0[1],len(v0[-1]))
	if len(v0[-1])>10:
		loc1.append(v0)
	else:
		for va in v0[-1]:
			vabord.append(va)
	print (len(loc1),'block num',len(vabord))
	
	n0=0
	dict3={}	
	for v3 in v2:
		if v3[1] not in vabord:
			dict3[v3[1]]=n0
			n0+=1
		#else:
		#	print (vabord[:10])
		#	print (v3[1])
	print (len(dict3),num-len(vabord))
	fig,ax=plt.subplots(figsize=(15,15))
	n=0
	heat1=[]	
	nlen=num+((len(loc1)-1)*20)-len(vabord)
	print (nlen)
	while n<=nlen:
		n1=[0]*nlen
		heat1.append(n1)
		n+=1
	T=1
	print (len(heat1),num+((len(loc1)-1)*20)-len(vabord))
	dictlo={}
	dlen=[]
	dnum=[]
	for k,v in dict21.items():
		dlen.append(len(k))
		dnum.append(len(set(v)))
		for j in v:
			if j.split('\t')[0] in vabord:
				#print ('contunh')
				#sys.exit()
				continue	
			x=dict3[j.split('\t')[0]]
			xa=0
			for lo in loc1:
				lon=str(lo[0])+"-"+str(lo[1])
				if j.split('\t')[0] in lo[-1]:
					lon=str(lo[0])+"-"+str(lo[1])
					if lon in dictlo:
						listlon=dictlo[lon]
					else:
						listlo=0
					listlo=len(set(lo[-1]))
					dictlo[lon]=listlo
					break
				else:
					xa+=20
			for j1 in v:
				if j1.split('\t')[0] in vabord:
					continue
				y=dict3[j1.split('\t')[0]]
				ya=0
				for lo in loc1:
					if j1.split('\t')[0] in lo[-1]:
						lon=str(lo[0])+"-"+str(lo[1])
						if lon in dictlo:
							listlo=dictlo[lon]
						else:
							listlo=0
						listlo=len(set(lo[-1]))
						dictlo[lon]=listlo
						break
					else:
						ya+=20
				#print (y+ya-1,x+xa-1)
				#if xa==40 and ya==40:
				#	print (j1.split('\t')[0],dict4[j1.split('\t')[0]])
				#print (dict4[j1.split('\t')[0]],dict4[j.split('\t')[0]],dict4[j1.split('\t')[0]]==dict4[j.split('\t')[0]],T)
				heat1[y+ya-1][x+xa-1]=T#int(dict4[j1.split('\t')[0]])
				#print (dict4[j1.split('\t')[0]])
		T+=1
	
	print (len(heat1))
	#cmap=mpl.cm.viridis
	#norm=mpl.colors.Normalize(vmax=6,vmin=-1)
	#cmap=mpl.Pastel1
	print (np.linspace(-3.5, 3.5, 8))
	print (dict(ticks=np.arange(-3,4)))
	print (np.linspace(-0.5,5.5,7))
	print (dict(ticks=np.arange(0,6)))
	#print (heat1.count(0.0))
	#cbar_kw=dict(ticks=np.arange(0,6))
	print (heat1[0][:5])
	#im,cbar = heatmap(heat1,ax=ax,cbarlabel="Type of CentGm-4",cmap=mpl.colors.ListedColormap(["#ffffff","#1f78b4","#33a02c","#e31a1c","#ff7f00","#b15928"],name=[0,1,2,3,4,5],N=6))#,norm=colors.NoNorm())
	im,cbar = heatmap(heat1,ax=ax,cbarlabel="Type of CentGm-4",cmap=mpl.colors.ListedColormap(["#ffffff","#d62728","#2ca02c","#ff7f0e","#1f77b4","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"],name=[0,1,2,3,4,5,6,7,8,9,10],N=11))
	#im,cbar = heatmap(heat1,ax=ax,cbarlabel="Type of CentGm-4",cmap="Pastel1_r")
	
	#im,cbar = heatmap(heat1,ax=ax,cmap="cubehelix_r",cbarlabel="Identical of CentGm-4")
	listlo=[]
	for k,v in dictlo.items():
		ks,ke=k.split('-')
		listlo.append([int(ks),int(ke),str(ks)+"-"+str(ke),v])
	listlo1=sorted(listlo,key=(lambda x:x[0]))
	xn=[]
	xnum=[]
	n1=0
	xlo=[]
	for v in listlo1:
		print (v)
		xnum.append(n1+int(v[-1]/2))
		xn.append(v[-2])
		print (v,n1+int(v[-1]/2))
		n1=n1+int(v[-1])+20
		xlo.append([int(v[-2].split('-')[0]),int(v[-2].split('-')[1])])
	print (xnum)
	print (xlo)
	ax.set_xticks(xnum,labels=xn)
	#ax.set_xticklabels(xnum)
	ax.set_yticklabels([])
	#fig.colorbar(im, ax=ax)
	figname=name+"-CentGm-4.pdf"
	fig.savefig(figname)
	plt.close()
	
	fig,ax=plt.subplots(2,1,sharex=True,figsize=(30.7,8.27))
	lables=['Length','Number']
	c=["#d62728","#2ca02c","#ff7f0e","#1f77b4","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"]
	x=range(2)
	width = 0.35 
	xn=[]
	for i in range(1,11,1):
		x=[i]
		y=dlen[i-1]
		clor=c[i-1]
		y1=dnum[i-1]
		xn.append("Type"+str(i))	
		rects1=ax[0].bar([i],[y],label="Type"+str(i),align="center",width=width,color=clor)
		ax[0].bar_label(rects1, padding=3)
		rects2=ax[1].bar([i],[y1],label="Type"+str(i),align="center",width=width,color=clor)
		ax[1].bar_label(rects2, padding=3)
	ax[0].set_xticks(range(1,11,1),labels=xn)
	
	ax[0].set_ylabel('Length of CentGm-4')
	ax[1].set_ylabel('Number of CentGm-4')
	figname="CentGm-4-"+name+".Type.pdf"
	fig.savefig(figname)
	plt.close()
	print ("welldone")
