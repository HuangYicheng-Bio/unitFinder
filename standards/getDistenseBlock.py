import re
import Bio
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import getopt
import sys
import numpy as np
import matplotlib
import matplotlib.ticker as ticker
def usage():
	print ("--unit *block.units.fa")
	print ("--TE *.RM.txt")
	print ("--name ID")
	print ("--sp sp name")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["name=","unit=","TE=","start=","end=",'chr=','all='])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq2=''
	TE=''
	start=''
	end=''
	name=''
	Chr=''
	all1=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--unit"):
			seq2=value
		elif opt in ("--name"):
			name=value
		elif opt in ("--start"):
			start=int(value)
		elif opt in ("--end"):
			end=int(value)
		elif opt in ("--TE"):
			TE=value
		elif opt in ("--chr"):
			Chr=value
		elif opt in ("--all"):
			all1=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	getUniqBlock(seq2,TE,name,start,end,Chr,all1)
def getDistenseBlock(pwd2,pwd3,name,start,end,Chr,all1,covpwd,pwdu1,pwdu4):
	dict1={}
	dict4={}
	if covpwd=='':
		#fig, axt = plt.subplots(4,1,sharex=True)
		#gs = axt[1].get_gridspec()
		#ax1=fig.add_subplot(gs[1:-1, :],sharex=axt[-1])
		#ax1.tick_params(axis="x", labelbottom=False)
		#axt[2].remove()
		#axt[1].remove()
		#ax=[axt[0],ax1,axt[-1]]
		fig, ax= plt.subplots(3,1,sharex=True)
		z=1
	else:
		#fig, axt = plt.subplots(5,1,sharex=True)
		#gs = axt[2].get_gridspec()
		#ax1=fig.add_subplot(gs[2:-1, :],sharex=axt[-1])
		#ax1.tick_params(axis="x", labelbottom=False)
		#axt[3].remove()
		#axt[2].remove()
		#ax=[axt[0],axt[1],ax1,axt[-1]]
		fig, ax= plt.subplots(4,1,sharex=True)
		z=0
	print (ax)
	dictseq={}
	max1=0
	for gseq in SeqIO.parse(all1,'fasta'):
		if gseq.seq in dictseq:
			listt=dictseq[gseq.seq]
		else:
			listt=[]
		listt.append(gseq.id)
		dictseq[gseq.seq]=listt
		if max1<len(listt):
			max1=len(listt)
	print (max1)
	dict1seq={}
	dictx={}
	x0=1
	for gseq in SeqIO.parse(pwdu1,'fasta'):
		if gseq.seq in dict1seq:
			listt=dict1seq[gseq.seq]
		else:
			listt=[]
			dictx[gseq.seq]=x0
			x0+=2
		listt.append(gseq.id)#.split(".unit-")[0])
		dict1seq[gseq.seq]=listt
	listt1=[]
	for k,v in dict1seq.items():
		listt1.append([k,len(v)])
	listt2=sorted(listt1,key=(lambda x:x[1]),reverse=True)
	dict1seqtop={}
	T=1
	list1seqtop=[]
	print ("CentGm-1")
	for i in listt2[:10]:
		dict1seqtop["Type"+str(T)]=dict1seq[i[0]]
		T+=1
		list1seqtop=list1seqtop+dict1seq[i[0]]
		print (len(dict1seq[i[0]]),i,len(i[0]))
	dict4seq={}
	
	for gseq in SeqIO.parse(pwdu4,'fasta'):
		if gseq.seq in dict4seq:
			listt=dict4seq[gseq.seq]
		else:
			listt=[]
			dictx[gseq.seq]=x0
			x0+=2
		listt.append(gseq.id)#.split(".unit-")[0])
		dict4seq[gseq.seq]=listt
	listt1=[]
	for k,v in dict4seq.items():
		listt1.append([k,len(v)])
	listt2=sorted(listt1,key=(lambda x:x[1]),reverse=True)
	dict4seqtop={}
	T=1
	list4seqtop=[]
	print ("CentGm-4")
	for i in listt2[:10]:
		dict4seqtop["Type"+str(T)]=dict4seq[i[0]]
		print (len(dict4seq[i[0]]),i,len(i[0]))
		T+=1
		list4seqtop=list4seqtop+dict4seq[i[0]]
	category_colors4 = plt.colormaps['twilight_shifted'](np.linspace(0, 1, max1))
	#category_colors4 = plt.colormaps['cividis'](range(0, max1+1, 1))
	print (len(category_colors4))
	
	si=start
	dicttype={}
	tx=[]
	#for i in range(start,end+10000,10000):
	#	if i>si:
	#		k=str(si)+"-"+str(i)
	#		dicttype[k]=[]
	#		tx.append([si,i])
	#	si=i
	#if i <end:
	#	k=str(si)+"-"+str(i)
	#	dicttype[k]=[]
	#	tx.append([si,i])
	locx=[]
	for gseq in SeqIO.parse(all1,'fasta'):
		ID=gseq.description.split('\t')[-2].split('-')
		x=[]
		y=[]
		x.append(int(ID[0])+start)
		y.append(len(gseq.seq))
		x.append(int(ID[1])+start)
		y.append(len(gseq.seq))
		locx.append([int(ID[0])+start,gseq.seq])
	#	for k,v in dicttype.items():
	#		s,e=k.split("-")
	#		st=int(ID[0])+start
	#		et=int(ID[1])+start
	#		xt=range(int(s),int(e)+1,1)
	#		yt=range(st,et+1,1)
	#		t=set(xt).intersection(set(yt))
	#		if len(t)!=0:
	#			listt2=dicttype[k]
	#			listt2.append(gseq.seq)
	#			dicttype[k]=listt2
		#if gseq.id.split(".unit-")[0] in list1seqtop:
		if gseq.id in list1seqtop:
			ax[3-z].axvspan(int(ID[0])+start,int(ID[1])+start, facecolor='green',alpha=0.5)
		#if gseq.id.split(".unit-")[0] in list4seqtop:
		if gseq.id in list4seqtop:
			ax[3-z].axvspan(int(ID[0])+start,int(ID[1])+start, facecolor='red',alpha=0.5)
		#ax[3-z].plot(x, y,color=matplotlib.colors.to_hex(category_colors4[len(dictseq[gseq.seq])/float(max1)]))
		ax[3-z].plot(x, y,lw=.1,color=matplotlib.colors.to_hex(category_colors4[len(dictseq[gseq.seq])-1]))
	#ax[3-z].set_title('CentGm Units', size="small")
	#ax[3-z].text(0.55, 0.1, "The top 10 CentGm-4 Units with the highest abundance",ha="right", va="top",bbox=dict(boxstyle="square",ec="red",fc="red",alpha=0.5))
	#ax[3-z].text(0.55, 0.2, "The top 10 CentGm-1 Units with the highest abundance",ha="right", va="top",bbox=dict(boxstyle="square",ec="green",fc="green",alpha=0.5))
	ax[3-z].plot([],[],color="green",label="The top 10 highest abundance CentGm-1 Units")
	ax[3-z].plot([],[],color="red",label="The top 10 highest abundance CentGm-4 Units")
	ax[3-z].legend(ncols=2,fontsize="xx-small",bbox_to_anchor=(0, 1),loc='lower left',labelcolor="linecolor")
	ax[3-z].set_ylabel("Unit Lenght (bp)")
	maxn=0
	#ax[3-z].colorbar(ax=ax[3-z], label='Number of Different Units')
	ty=[]
	tx0=[]
	
	#for k in tx:
	#	s0,e0=k
	#	k1=str(s0)+"-"+str(e0)
	#	tx0.append(s0)
	#	ty.append(len(set((dicttype[k1]))))
	#	tx0.append(e0)
	#	ty.append(len(set((dicttype[k1]))))
	locx1=sorted(locx,key=(lambda x:x[0]))
	locxt=[]
	lx=[]
	for i in locx1:
		if len(locxt)<50:
			locxt.append(i[1])
			lx.append(i[0])
		else:
			ty.append(len(set(locxt)))
			tx0.append((max(lx)+min(lx))/2.0)
			locxt=[]
			lx=[]
			locxt.append(i[1])
			lx.append(i[0])
	if len(lx)!=0:
		ty.append(len(set(locxt)))
		tx0.append(np.mean(lx))
	par1 = ax[0].twinx()
	par1.plot(tx0, ty,color="black",lw=0.5,label="Type number per 50 Units")
	par1.legend(fontsize="xx-small",loc="upper right",labelcolor="linecolor")
	par1.set_ylabel("Type number")
	par1.yaxis.get_label().set_color("black")
	for gseq in SeqIO.parse(pwd2,'fasta'):
		ID=gseq.description.split('\t')[-1]
		if ID in dict1:
			list1t=dict1[ID]
		else:
			list1t=[]
		list1t.append(gseq)
		dict1[ID]=list1t
	for k,v in dict1.items():
		dict2={}
		for gseq in v:
			if gseq.seq in dict2:
				list1t=dict2[gseq.seq]
			else:
				list1t=[]
			loc=gseq.description.split('\t')[-3].split('-')[0]
			list1t.append(int(loc))
			dict2[gseq.seq]=list1t
		len1=[]
		t=[]
		lenp=[]
		for k1,v1 in dict2.items():
			if len(v1)>5: #cp number
				v2=sorted(v1)
				s0=v2[0]
				for i1 in v2[1:]:
					a=i1-s0
					len1.append([a,i1,s0])
					s0=i1
		if len(len1)!=0:
			list2=sorted(len1,key=(lambda x:x[1]))
			if len(list2)>maxn:
				maxn=len(v)
	dict3={}
	lenall=[]
	list2=[]
	#category_colors4 = plt.colormaps['Set3'](np.linspace(0, len(dict1)-1, len(dict1)))
	norm = mpl.colors.Normalize(vmin=0, vmax=len(dict1))
	print (np.linspace(0, len(dict1), len(dict1)))
	c0=0
	boxt=[]
	boxc=[]
	boxp=[]
	par2 = ax[2-z].twinx()
	for k,v in dict1.items():
		dict2={}
		colorR=plt.colormaps['Set3'](c0)
		
		#print (colorR)

		for gseq in v:
			if gseq.seq in dict2:
				list1t=dict2[gseq.seq]
			else:
				list1t=[]
			loc=gseq.description.split('\t')[-3].split('-')[0]
			list1t.append(int(loc))
			dict2[gseq.seq]=list1t
		len1=[]
		t=[]
		lenp=[]
		x=1
		xd=2
		for k1,v1 in dict2.items():
			al=[]
			lent=[]
			if len(v1)>5: #cp number
				v2=sorted(v1)
				s0=v2[0]
				for i1 in v2[1:]:
					a=i1-s0
					#if len(k1)*100>a:
					al.append(a)
					len1.append([a,s0,i1])
					lent.append([a,s0,i1])
					lenall.append([a,s0,i1])
					s0=i1
			if len(lent)!=0:
				list2=sorted(lent,key=(lambda x:x[1]))
				t=[]
				lenp=[]
				for i in list2:
					#x=dictx[k1]
					t.append(i[1]+start)
					#lenp.append(i[0])
					lenp.append(x)
					t.append(i[2]+start)	
					#lenp.append(i[0])
					lenp.append(x)
				ax[2-z].plot(t, lenp,linestyle='dashed',color=matplotlib.colors.to_hex(colorR),linewidth=0.1)
				ax[2-z].scatter(t,lenp,color=matplotlib.colors.to_hex(colorR),s=0.1)
				x+=xd
				

		if len(len1)!=0:
			dict3[k]=len1
			list2=sorted(len1,key=(lambda x:x[1]))
			t0=[]
			lociD=[]
			list2t=[]
			tl0=[]
			lociDl=[]
			for i in list2:
				t0.append(i[0])
				lociD.append(i[1]+start)
				lociD.append(i[2]+start)
				tl0.append(i[0])
				lociDl.append(i[1]+start)
		
				list2t.append([i[0],i[1]])
				list2t.append([i[0],i[1]])
		

			print (min(lociD),0,max(lociD)-min(lociD),max(t0)-0,gseq.description,len(v))

			wth=(end-start)/50
			bplot=par2.boxplot(t0,showfliers=False,patch_artist=True,medianprops={"color": "red", "linewidth": 0.5},boxprops={"facecolor": colorR, "edgecolor":"black","linewidth": 0.5},whiskerprops={"color": 'black', "linewidth": 0.5},capprops={"color": 'black', "linewidth": 0.5},widths=wth,positions=np.array([(min(lociD)+max(lociD))/2.0]))
			boxt.append(t0)
			boxc.append(colorR)
			par2.set_ylabel("Distanse (bp)")
			boxp.append((min(lociD)+max(lociD))/2.0)
			print (list(range(1,2)))
			#t=[]
			#lenp=[]
			#t0=[]
			#lenp0=[]
			#list2t1=sorted(list2t,key=(lambda x:x[1]))
			#for i in list2:
			#	t=[]
			#	lenp=[]
			#	t.append(i[1]+start)
			#	lenp.append(i[0])
			#	t.append(i[2]+start)
			#	lenp.append(i[0])
			#	t0.append(i[1]+start)
			#	lenp0.append(i[0])
			#	t0.append(i[2]+start)
			#	lenp0.append(i[0])
			#	ax[2-z].plot(t, lenp,linestyle='dashed',color=matplotlib.colors.to_hex(colorR),alpha=float(len(list2)/maxn),linewidth=0.5)
			#ax[2-z].scatter(t0,lenp0,color=matplotlib.colors.to_hex(colorR),s=0.5)
			c0+=1
			if c0>=11:
				c0=0
	ax[2-z].set_title('Distanse within Blocks', size="small")
	ax[2-z].set_ylabel("Units")
	file1=open(pwd3,'r')
	list1=[]
	t=[]
	tall=[]
	for row in file1:
		row1=row.rstrip().split(' ')
		row2=[]
		for r in row1:
			if r!='':
				row2.append(r)
		tall.append(int(row2[6])-int(row2[5]))
	mean1=np.mean(tall)
	std1=np.std(tall)
	t3=3*std1
	file1.seek(0,0)
	list1in=[]
	tin=[]
	for row in file1:
		row1=row.rstrip().split(' ')
		if 'Simple_repeat' in row or 'Low_complexity' in row:
			continue
		row2=[]
		for r in row1:
			if r!='':
				row2.append(r)
		t0=int(row2[6])-int(row2[5])
		#if np.abs(t0-mean1)<=t3:
		if True:
			if "(0)" in row2:
				t.append(int(row2[5])-0.5)
				list1.append(0)
				t.append(int(row2[5]))
				list1.append(int(row2[6])-int(row2[5]))
				t.append(int(row2[6]))
				list1.append(int(row2[6])-int(row2[5]))
				t.append(int(row2[6])+0.5)
				list1.append(0)
			else:
				tin.append(int(row2[5])-0.5)
				list1in.append(0)
				tin.append(int(row2[5]))
				list1in.append(int(row2[6])-int(row2[5]))
				tin.append(int(row2[6]))
				list1in.append(int(row2[6])-int(row2[5]))
				tin.append(int(row2[6])+0.5)
				list1in.append(0)
	ax[1-z].plot(tin, list1in,lw=0.5,label="incomplete TE")
	ax[1-z].fill_between(tin,0,list1in, alpha=0.5)
	ax[1-z].plot(t, list1,lw=0.5,label="complete TE")
	ax[1-z].fill_between(t,0,list1, alpha=0.5)
	ax[1-z].set_title('TE location', size="small")
	ax[1-z].set_ylabel("TE Length (bp)")
	ax[1-z].legend(fontsize="xx-small",loc="upper right")
	figname1=name+".pdf"
	
	if covpwd!='':
		counts={}
		file2=open(covpwd,'r')
		for row in file2:
			row=row.rstrip()
			filet=open(row,'r')
			t=[]
			d=[]
			d0=[]
			for row1 in filet:
				row2=row1.rstrip().split('\t')
				if Chr in row2[0]:
					s=int(row2[1])
					e=int(row2[2])
					if s<=end and e>=start:
						t.append(s)
						d.append(int(row2[3]))
						t.append(e)
						d.append(int(row2[3]))
					#if int(row2[1])>=start and int(row2[1])<=end:
					#t.append(int(row2[1]))
					#d.append(int(row2[2]))
					#if int(row2[2])!=0:
		#			d0=d0+list(range(s,e+1,1))
		#d02=[]
		#for i in range(start,end,1):
		#	if i not in d0:
		#		d02.append(i)
		#if len(d02)!=0:
		#	d01=sorted(d02)
		#	d0t=[d01[0]]
		#	s0=d01[0]
		#	for di in d01[1:]:
		#		if di<=s0+1:
		#			d0t.append(di)
		#		else:
					#print (di<=s0+10,di)
		#			if len(d0t)>10000:
		#				print ("zero",min(d0t),max(d0t),(max(d0t)-min(d0t)))
		#				ax[0].axvspan(min(d0t),max(d0t), facecolor='grey',alpha=0.5)
		#				d0t=[di]
		#		s0=di
		#	print (len(d0))
			RNAseqID=row.split("\\")[-1].split(".")[0]
			ax[0].plot(t, d,lw=0.5,alpha=0.5,label='ChIP-seq '+RNAseqID)
			ax[0].fill_between(t,0,d, alpha=0.5)
			filet.close()
		#	RNAseqID=row.split("\\")[-1].split(".")[0]
		#	counts['RNA-seq '+RNAseqID]=[np.array(d),np.array(t)]
		file2.close()
		#ax[0].legend(fontsize="xx-small",loc="upper left",labelcolor="linecolor")
		ax[0].set_ylabel("Depth")
		ax[0].set_title('ChIP-seq and Type number', size="small")
		#bottom = np.zeros(len(d))
		#for k,v in counts.items():
		#	p=ax[0].bar(v[0],v[1],100,label=k, bottom=bottom)
		#	bottom+=v
		#ax[0].set_title('CENH3')
		#plt.title(name)
	#print (np.arange(start, end, 20))
	ax[0].set_xticks(np.arange(start, end, int((end-start)/5)))
	ax[0].xaxis.set_major_formatter(lambda x, pos: str(round(x/1000000,2))+" Mbp")
#	ax[0].xaxis.set_major_locator(ticker.MultipleLocator(1000000))
#	ax[0].xaxis.set_minor_locator(ticker.MultipleLocator(100000))
	#ax[0].xaxis.set_major_formatter(lambda x, pos: str(round(x/1000000))+" Mbp")
	ax[3-z].set_ylabel(Chr)
	fig.tight_layout()
	fig.savefig(figname1)
	plt.close()
#getoptions()
