import re
import os
import Bio
from Bio import SeqIO
import getopt
import sys
def usage():
	print ("--in *.fa")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name="])
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
		elif opt in ("--name"):
			name=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,name

def getunit(roundNumPwd,name,roundNum):
	dictcosname={}
	p=name+".units.fa"
	fileo=open(p,'w')
	commandline='cat '+roundNumPwd+"/*.txt.html >temp"
	os.system(commandline)
	file1t=open('temp','r')
	dict1={}
	ID1=''
	jr=0
	ID1=''
	unitdicts={}
	seqlist=[]
	for row in file1t:
		if 'Sequence:' in row:
			row1=row.rstrip().split(' ')
			r=[]
			for i in row1:
				if i!='':
					r.append(i)
			ID=r[1].split('\t')[0]
			if ID not in seqlist:
				seqlist.append(ID)
				jr=0
		elif '<A NAME="' in row:
			jr+=1
			ID1=ID+'_Repeat_'+str(roundNum)+'-'+str(jr)
			dict1[ID1]={}
			seq1=[]
			ca=row.split('">')[0].split('<A NAME="')[1]
			ca1=ca.split(',')
			ca2=ID+"-"+','.join(ca1[:2])
			if ca2 in dictcosname:
				print (ca2)
				print (ID1,dictcosname[ca2])
				sys.exit()
			else:
				dictcosname[ca2]=ID1
		elif row[0]==' ' or row=='\n':
			if ID1!='':
				row1=row.rstrip().split(' ')
				r=[]
				for i in row1:
					if i!='' and '*' not in row:
						r.append(i)
				if 'Indices' in row:
					r=[]
				if 'Period' in row:
					len1=int(r[4].split('.')[0])
					clen=int(r[-1])
					r=[]
				if len(r)>2:
					r1='_'.join(r[1:])
					r0=[r[0],r1]
					seq1.append(r0)
				else:
					seq1.append(r)
		elif 'Statistics' in row:
			j=0
			seq21=[]
			seq22=[]
			j0=0
			tag=1
			while j0<len(seq1):
				i=seq1[j0]
				if i!=[]:
					if j==0:
						if j0-1>=0 and j0+2<len(seq1):
							if seq1[j0-1]==[] and seq1[j0+2]==[] and seq1[j0+1]!=[]:
								seq21.append(i)
								j=1
								tag=1
					else:
						if j0+1<len(seq1) and j0-2>=0:
							if seq1[j0+1]==[] and seq1[j0-2]==[] and seq1[j0-1]!=[]:
								seq22.append(i)
								j=0
								tag=0
				j0+=1
			if len(seq21)!=len(seq22):
				print (len(seq21),len(seq22))
				print (seq22[0],seq22[1])
				print (seq21[0],seq21[1])
				print ('error')
				sys.exit()
			elif len(seq21)!=0:
				i=0
				seq1t=''
				seq2t=''
				seq31=[]
				seq32=[]
				t=1
				t0=[]
				t1=[]
				while i<len(seq21):
					if seq1t=='':
						seq1t=seq1t+seq21[i][1]
						seq2t=seq2t+seq22[i][1]
						t=int(seq22[i][0])
						t0.append(int(seq21[i][0]))
						t1.append(int(seq22[i][0]))
					else:
						if t<int(seq22[i][0]):
							seq1t=seq1t+seq21[i][1]
							seq2t=seq2t+seq22[i][1]
							t=int(seq22[i][0])
							t0.append(int(seq21[i][0]))
							t1.append(int(seq22[i][0]))
						else:
							seq31.append([seq1t,t0])
							seq32.append([seq2t,t1])
							seq1t=''
							seq2t=''
							seq1t=seq1t+seq21[i][1]
							seq2t=seq2t+seq22[i][1]
							t0=[]
							t1=[]
							t0.append(int(seq21[i][0]))
							t1.append(int(seq22[i][0]))
							t=int(seq22[i][0])
					i+=1
				seq31.append([seq1t,t0])
				seq32.append([seq2t,t1])
				i=0
				if len(seq32)-len1==0:
					seq42=seq32
					seq41=seq31
				else:
					i=len1-len(seq32)
					seq42=seq32[:i]
					seq41=seq31[:i]
				ju=1
				i=0
				while i<len(seq41):
					seq1o,loc=seq41[i]
					seq2o,loc2=seq42[i]
					if '_' not in seq1o or '_' not in seq2o:
						seqo2=seq2o.replace('-','')
						seqo=seq1o.replace('-','')
						if len(seq2o)<clen:
							lo='>'+ID1+'.unit-'+str(ju)+"\t"+str(min(loc))+"-"+str(min(loc)+len(seqo))+"\t"+str(min(loc2))+"-"+str(min(loc2)+len(seqo2))+"\n"+seqo+"\n"
						else:
							lo='>'+ID1+'.unit-'+str(ju)+"\t"+str(min(loc))+"-"+str(min(loc)+len(seqo))+"\t"+str(min(loc2))+"-"+str(min(loc2)+clen)+"\n"+seqo+"\n"
						if ID1 in unitdicts:
							listt=unitdicts[ID1]
						else:
							listt=[]
						listt.append(ID1+'.unit-'+str(ju))
						unitdicts[ID1]=listt
						fileo.writelines(lo)
						ju+=1
						
					else:
						seqo1=seq1o.split('_')
						seqo2=seq2o.split('_')
						j0=0
						jlen=0
						while j0<len(seqo1):
							seqo11=seqo1[j0].replace('-','')
							seqo21=seqo2[j0].replace('-','')
							if len(seqo2)<clen:
								lo='>'+ID1+'.unit-'+str(ju)+"\t"+str(min(loc)+jlen)+"-"+str(min(loc)+jlen+len(seqo11))+"\t"+str(min(loc2))+"-"+str(min(loc2)+len(seqo21))+"\n"+seqo11+"\n"
							else:
								lo='>'+ID1+'.unit-'+str(ju)+"\t"+str(min(loc)+jlen)+"-"+str(min(loc)+jlen+len(seqo11))+"\t"+str(min(loc2))+"-"+str(min(loc2)+clen)+"\n"+seqo11+"\n"
							if ID1 in unitdicts:
								listt=unitdicts[ID1]
							else:
								listt=[]
							listt.append(ID1+'.unit-'+str(ju))
							unitdicts[ID1]=listt
							j0+=1
							ju+=1
							jlen=jlen+len(seqo11)
							fileo.writelines(lo)

					i+=1
	file1t.close()
	fileo.close()
	commandline='rm temp'
	os.system(commandline) 
	return p,unitdicts,dictcosname

def getdot(pwdtrf,roundNum,name,dlocpro,dictcosname):
	file1=open(pwdtrf,'r')
	pwd1dot=name+'.consensus.fa'
	pwd2dot=name+'.all.fa'
	file2=open(pwd1dot,'w')
	file3=open(pwd2dot,'w')
	j=1
	repdict={}
	SeqID=''
	seqlist=[]
	for row in file1:
		if 'Sequence:' in row:
			if SeqID!='':
				repdict[SeqID]=replist
			row1=row.rstrip().split(' ')
			SeqID=row1[-1].split('\t')[0]
			replist=[]
			if SeqID not in seqlist:
				seqlist.append(SeqID)
				j=1
		row1=row.rstrip().split(' ')
		if len(row1)>=12:
			s0=dlocpro[SeqID]-1
			ca1=SeqID+"-"+row1[0]+"--"+row1[1]+","+row1[2]
			IDc=dictcosname[ca1]
			row2=[str(int(row1[0])+s0),str(int(row1[1])+s0)]
			row1=row2+row1
			l='>'+IDc+"\t"+str(len(row1[-2]))+"-"+row1[5]+"\t"+"\t".join(row1[:-2])+"\n"+row1[-2]+"\n"
			file2.writelines(l)
			l='>'+IDc+"\t"+str(len(row1[-2]))+"-"+row1[5]+"\t"+"\t".join(row1[:-2])+"\n"+row1[-1]+"\n"
			file3.writelines(l)
			replist.append(IDc)
			j+=1
	if SeqID!='':
		repdict[SeqID]=replist

	file1.close()
	file2.close()
	file3.close()

	return pwd1dot,pwd2dot,repdict

def checkrep(repdict,pwdpro,pwdprounit,conspwd,fileo1,fileo2):
	vflag=[]
	for k,v in repdict.items():
		if len(v)==0:
			for gseq in SeqIO.parse(pwdpro,'fasta'):
				if gseq.id==k:
					l='>'+gseq.description+"\n"+gseq.seq+"\n"
					fileo1.writelines(l)
					for gseq1 in SeqIO.parse(pwdprounit,'fasta'):
						k1=gseq1.id.split('.unit-')[0]
						if k1==gseq.id:
							l='>'+gseq1.description+"\n"+gseq1.seq+"\n"
							fileo2.writelines(l)
		else:
			vflag+=v
	return vflag

def nextinput(vflag,conspwd,unitspwd,unitdict,name):
	pwd1=name+'.next.fa'
	pwd2=name+'.next.unit.fa'
	
	file1=open(pwd1,'w')
	for gseq in SeqIO.parse(conspwd,'fasta'):
		if gseq.id in vflag:
			l='>'+gseq.description+"\n"+gseq.seq+"\n"
			file1.writelines(l)
	v=[]
	for i in vflag:
		v=v+list(unitdict[i])
	file2=open(pwd2,'w')
	for gseq in SeqIO.parse(unitspwd,'fasta'):
		vid=gseq.id.split('.unit-')[0]
		if vid in vflag:
			l='>'+gseq.description+"\n"+gseq.seq+"\n"
			file2.writelines(l)
	file1.close()
	file2.close()
	return pwd1,pwd2
def unified(conspwd,unitspwd,tag,roundNum):
	file1=open('temp','w')
	d1={}
	for gseq in SeqIO.parse(conspwd,'fasta'):
		row=gseq.description.split('\t')
		if tag==0:
			s=int(row[2])
		else:
			s=int(row[-1].split('-')[0])
		d1[gseq.id]=s
	for gseq in SeqIO.parse(unitspwd,'fasta'):
		if tag==1:
			ID0=gseq.id.split('_Repeat_')
			ID01='_Repeat_'.join(ID0[:-1])
			ID=ID01
		else:
			ID01=gseq.id
			ID=ID01.split('.unit')[0]
		s=int(gseq.description.split('\t')[1].split('-')[0])
		e=int(gseq.description.split('\t')[1].split('-')[1])
		if roundNum!=1:
			s0=d1[ID]-1
		else:
			s0=0
		l='>'+gseq.description+"\t"+str(s0+s)+"-"+str(s0+e)+"\n"+gseq.seq+"\n"
		file1.writelines(l)
	file1.close()
	commondline='mv temp '+unitspwd
	os.system(commondline)
	return d1
def changeunitID(unitspwd,unitdictpro):
	unitdict1={}
	file1=open('temp','w')
	for gseq in SeqIO.parse(unitspwd,'fasta'):
		k=gseq.id.split('_Repeat_')
		k1='_Repeat_'.join(k[:-1])
		
		u1=k[-1].split('.unit-')[-1]
		u2=k[-1].split('.unit-')[-2].split('-')[-1]
		i=0
		ID=k1+'-'+u2+","+u1
		ID2=gseq.description.split('\t')
		l='>'+ID+"\t"+"\t".join(ID2[1:])+"\n"+gseq.seq+"\n"
		file1.writelines(l)
	file1.close()

	for k0,v0 in unitdictpro.items():
		kt=k0.split('_Repeat_')
		kt1='_Repeat_'.join(kt[:-1])
		k1=kt1.split('.unit-')[0]
		lt1=[]
		for k00 in v0:
			k01=k00.split('_Repeat_')
			k011='_Repeat_'.join(k01[:-1])
			u1=k01[-1].split('.unit-')[-1]
			ID=k011+'-'+u1
			lt1.append(ID)
		unitdict1[k1]=lt1
	commondline='mv temp '+unitspwd
	os.system(commondline)
	return unitdict1
	

def trf(pwdpro,pwdprounit,dlocpro,roundNum,name,fileo1,fileo2):
	nametrf=pwdpro.split('/')[-1].split('.fa')[0]
	pwdtrf=nametrf+'.fa.2.7.7.80.10.50.500.dat'
	print (pwdtrf)
	if os.path.exists(pwdtrf)!=True or os.path.getsize(pwdtrf)==0:
		commandline='trf '+pwdpro+' 2 7 7 80 10 50 500 -f -d -m'
		os.system(commandline)
		commandline='mkdir round'+str(roundNum)
		os.system(commandline)
		commandline='mv *.html ./round'+str(roundNum)
		os.system(commandline)
	roundpwd='round'+str(roundNum)
	unitspwd,unitdict,dictcosname=getunit(roundpwd,name,roundNum)

	conspwd,allpwd,repdict=getdot(pwdtrf,roundNum,name,dlocpro,dictcosname)
	dloc=unified(conspwd,unitspwd,0,roundNum)
	
	vflag=checkrep(repdict,pwdpro,pwdprounit,conspwd,fileo1,fileo2)
	
	unitspwdpro=''
	if pwdprounit!='':
		pwdtunit=name+".pro.units.fa"
		filet=open(pwdtunit,'w')
		for gseq in SeqIO.parse(conspwd,'fasta'):
			ID1=gseq.id.split('_Repeat_')
			ID2='_Repeat_'.join(ID1[:-1])
			for gseq1 in SeqIO.parse(pwdprounit,'fasta'):
				id1=gseq1.id.split('.unit-')
				a=gseq1.description.split('\t')
				if ID2==id1[0]:
					u='.unit-'.join(id1[1:])
					l='>'+gseq.id+'.unit-'+u+"\t"+"\t".join(a[1:])+"\n"+gseq1.seq+"\n"
					filet.writelines(l)
		filet.close()
		nametrf=pwdtunit.split('/')[-1].split('.fa')[0]
		pwdtrf=nametrf+'.fa.2.7.7.80.10.50.500.dat'
		if os.path.exists(pwdtrf)!=True or os.path.getsize(pwdtrf)==0:
			commandline='trf '+pwdtunit+' 2 7 7 80 10 50 500 -f -d -m'
			os.system(commandline)
			commandline='mkdir unit_round'+str(roundNum)
			os.system(commandline)
			commandline='mv *.html ./unit_round'+str(roundNum)
			os.system(commandline)
		roundpwd='unit_round'+str(roundNum)
		name1=name+".pro.units"
		unitspwdpro,unitdictpro,dictcosname=getunit(roundpwd,name1,roundNum)
		dloc1=unified(pwdtunit,unitspwdpro,1,roundNum)
		unitdictpro=changeunitID(unitspwdpro,unitdictpro)
	#print (vflag[:5],roundNum)
	if unitspwdpro=='':
		nextpwd,nextunitpwd=nextinput(vflag,conspwd,unitspwd,unitdict,name)
		return unitspwd,conspwd,allpwd,vflag,nextpwd,nextunitpwd,dloc
	else:
		for k,v in unitdict.items():
			if k in unitdictpro:
				lv=unitdictpro[k]
			else:
				lv=[]
			lv=lv+v
			unitdictpro[k]=lv
		pwddunit=name+".pro.all.units.fa"
		commandline='cat '+unitspwd+" "+unitspwdpro+" > "+pwddunit
		os.system(commandline)
		nextpwd,nextunitpwd=nextinput(vflag,conspwd,pwddunit,unitdictpro,name)
		return pwddunit,conspwd,allpwd,vflag,nextpwd,nextunitpwd,dloc

def secondary(pwd1,name):

	pwdo1=name+".Consensus.fa"
	pwdo2=name+".Units.fa"
	fileo1=open(pwdo1,'w')
	fileo2=open(pwdo2,'w')
	name=pwd1.split('/')[-1].split('.fa')[0]
	pwdprounit=''
	pwdpro=pwd1
	dproloc={}
	for gseq in SeqIO.parse(pwd1,'fasta'):
		dproloc[gseq.id]=1
	roundNum=1
	flag=1
	while flag==1:
		name1=name+'-'+str(roundNum)
		unitspwd,conspwd,allpwd,vflag,nextpwd,nextunitpwd,dloc=trf(pwdpro,pwdprounit,dproloc,roundNum,name1,fileo1,fileo2)
		if vflag==[]:
			flag=0
		else:
			pwdpro=nextpwd
			pwdprounit=nextunitpwd
			dproloc=dloc
		roundNum+=1
	fileo1.close()
	fileo2.close()
	print ("welldone")
