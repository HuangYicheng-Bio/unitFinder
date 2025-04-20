import re
import sys
import os
import Bio
from Bio import SeqIO
from mummer_2 import mummer2

#litt=['Chr04','Chr04','Chr04','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09','Chr10','Chr11','Chr12','Chr13','Chr14','Chr15','Chr16','Chr17','Chr18','Chr19','Chr20']
#litt=['Chr04']
#def completeGm1(Chrn,pwdt,pwd3,pwda,pwdc):
#for Chrn in litt:
def completeGm1(Chrn,units,type1,seq1,alnpwd):
	pwd2='./Standards/CentGm-1-'+Chrn+'.complete.txt'
	pwd4='./Standards/CentGm-1-'+Chrn+'.complete.fa'
	pwd5='./Standards/CentGm-1-'+Chrn+'.complete.next.txt'
	pwd6='./Standards/CentGm-1-'+Chrn+'.complete.coords'
	#pwd3='CentGm-1.typic.fa'
	file1=open(alnpwd,'r')
	dictu={}
	for gseq in SeqIO.parse(units,'fasta'):
		dictu[gseq.id]=gseq
	for gseq in SeqIO.parse(seq1,'fasta'):
		gseqall=gseq
	dict1={}
	dicta={}
	for row in file1:
		row1=row.rstrip().split('\t')
		if float(row1[6])<80:
			continue
		ID=row1[12]
		a=int(row1[2])
		b=int(row1[3])
		len1=int(row1[4])
		if a<b:
			f0='+'
		else:
			f0='-'
		if ID in dicta:
			aln,f=dicta[ID]
			if int(aln[4])<len1:
				aln=row1
				f=f0
				dicta[ID]=[aln,f]
		else:
			aln=row1
			f=f0
			dicta[ID]=[aln,f]
	file1.close()
	dicte={}
	dictue={}
	file2=open(pwd2,'w')
	file3=open(pwd4,'w')
	file4=open(pwd5,'w')
	filec=open(pwd6,'w')
	n0=0
	for k,v in dicta.items():
		gseq=dictu[k]
		print (k,n0)
		n0+=1
		loc=gseq.description.split('\t')[-1].split('-')
		aln,f=v
		print (loc,int(loc[0]),int(loc[1]))
		gseqt=gseqall.seq[int(loc[0])-1:int(loc[1])-1]
		print (gseqt)
		if gseqt!=gseq.seq:
			print ('seq error')
			print (aln)
			print (gseqt)
			print (gseq.seq)
			#sys.exit()
		s2=int(aln[2])
		e2=int(aln[3])
		print (s2,e2)
		sl=int(loc[0])-1
		el=int(loc[1])-1
		if f=='+':
			s=sl+s2-1
			e=s+(e2-s2+1)
		else:
			s=sl+e2-1
			e=s+(s2-e2+1)
	#	print (s,e,f)
	#	print (gseq.seq)
		#print (aln)
		if float(aln[9])<=95 or float(aln[10])<=95:
			s1=int(aln[0])-1+1
			e1=int(aln[7])-int(aln[1])+1
			if f=='+':
				sn=s-s1+1
				en=e+e1+1
				lentm=en-sn+1
				if lentm>len(gseq.seq):
					en=sn+len(gseq.seq)
			else:
				sn=s-e1+1
				en=e+s1-1
				lentm=en-sn+1
				if lentm>len(gseq.seq):
					sn=en-len(gseq.seq)
			#lentm=en-sn+1
			#if lentm>len(gseq.seq):
			#	en=sn+len(gseq.seq)
			if len(gseqall.seq[sn:en])>len(gseq.seq):
				print (len(gseqall.seq[sn:en]),len(gseq.seq))
				print (gseqall.seq[sn:en])
				print (gseq.seq)
				print ('lenerr')
				sys.exit()
		#	print (aln)
			ft=open('temp2.fa','w')
			l='>'+k+"\n"+gseqall.seq[sn:en]+"\n"
			print (gseqall.seq[sn:en])
			print ('After')
		#	print (l)
			ft.writelines(l)
			ft.close()
			#commandline='python ../mummer_1.py --seq1 '+pwd3+' --seq2 temp2.fa --name temp3'
			#os.system(commandline)
			mummer2(type1,'temp2.fa','temp3')
			ft=open('temp3.delta.filter.coords','r')
			alnt=[]
			Idg=[]
			match1=[]
			match2=[]
			for r in ft:
				print (r)
				r1=r.rstrip().split('\t')
				if float(r1[6])<60:
					continue
				IDr=r1[12]
				ar=int(r1[2])
				br=int(r1[3])
				len1r=int(r1[4])
				Idg.append(row1[6])
				match1+=range(int(r1[0])-1,int(r1[1])+1,1)
				if int(r1[2])<int(r1[3]):
					match2+=range(int(r1[2])-1,int(r1[3])+1,1)
				else:
					match2+=range(int(r1[3])-1,int(r1[2])+1,1)
				if ar<br:
					f0t='+'
				else:
					f0t='-'
				if alnt==[]:
					alnt=r1
					flt=f0t
				else:
					if int(alnt[4])<len1r:
						alnt=r1
						flt=f0t
			#print (len(set(match1))>=80,len(set(match2))>=80)
			#print (len(set(match1)),len(set(match2)),float(alnt[9]),alnt)
			#if float(alnt[4])>=80 and float(aln[5])>=80 and float(len(set(match1)))>=80 and float(aln[9])>=80:
			if len(set(match1))>=70 and len(set(match2))>=70:#and float(alnt[6])>=80:
				print ('Good')
				dicte[k]=[alnt,flt]
				dictue[k]=[gseq,sn,en]
				lo=str(sn)+"\t"+str(en)+"\t"+k+"\t"+str(alnt[0])+"\t"+str(alnt[1])+"\t-\t"+str(alnt[8])+"\n"
				file2.writelines(lo)
				if f=='+':
					lo='>'+k+"\t"+str(sn)+"-"+str(en)+"-"+str(alnt[0])+"-"+str(alnt[1])+"\t"+flt+"\n"+gseqall.seq[sn:en]+"\n"
				else:
					lo='>'+k+"\t"+str(sn)+"-"+str(en)+"-"+str(alnt[0])+"-"+str(alnt[1])+"\t"+flt+"\n"+gseqall.seq[sn:en].reverse_complement()+"\n"
				file3.writelines(lo)
				lo=str(sn)+"\t"+str(en)+"\t"+k+"\t"+r1[11]+"\t"+flt+"\t"+str(len(gseqall.seq[sn:en]))+"\t"+";".join(Idg)+"\t"+str(len(set(match1))/float(r1[7]))+"\n"
				file4.writelines(lo)
				filec.writelines("\t".join(alnt)+"\n")
			#else:
				#if 'GWHBWDJ00000001.1-15196401-16725249_Repeat_1-233' in k:
				#print (len(set(match1)),len(set(match2)),float(alnt[9]),alnt)
				#	print (f)
				#	print (gseqall.seq[sn:en])
				#	print ('>'+k+"\t"+str(sn)+"-"+str(en)+"-"+str(alnt[0])+"-"+str(alnt[1])+"\t"+flt+"\n"+gseqall.seq[sn:en].reverse_complement())
				#sys.exit()
			ft.close()
			commandline='rm temp*'
			os.system(commandline)
			#if f=="-":
				#print (gseqall.seq[sn:en].reverse_complement())
				#for gseq2 in SeqIO.parse(pwd3,'fasta'):
				#	print (gseq2.seq)
			
		else:
			dicte[k]=v
			dictue[k]=[gseq,sl,el]
			print (v)
			lo=str(sl)+"\t"+str(el)+"\t"+k+"\t"+str(aln[0])+"\t"+str(aln[1])+"\t-\t"+str(aln[8])+"\n"
			file2.writelines(lo)
			if f=='+':
				lo='>'+k+"\t"+str(sl)+"-"+str(el)+"-"+str(aln[0])+"-"+str(aln[1])+"\t"+f+"\n"+gseqall.seq[sl:el]+"\n"
			else:
				lo='>'+k+"\t"+str(sl)+"-"+str(el)+"-"+str(aln[0])+"-"+str(aln[1])+"\t"+f+"\n"+gseqall.seq[sl:el].reverse_complement()+"\n"
			file3.writelines(lo)
			lo=str(sl)+"\t"+str(el)+"\t"+k+"\t"+aln[11]+"\t"+f+"\t"+str(len(gseqall.seq[sl:el]))+"\t"+aln[6]+"\t"+str(float(aln[9])/100.0)+"\n"
			file4.writelines(lo)
			filec.writelines("\t".join(aln)+"\n")
	print ('Done')
	file2.close()
	file3.close()
	file4.close()
	file1.close()
	filec.close()
