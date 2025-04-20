import re
import os
import sys
import Bio
from Bio import SeqIO
pwd1='./ZH13.Type-TypeGroup/Type-7-typic.fasta'
pwd11='CentGm-4.fa'
#pwd11='../ref/Chr01-TypeGroup/Type-16-all.fasta'
pwd2='CentGm-4-stan.fa'
pwd3='CentGm-4-s-odd.fa'
file1=open(pwd2,'w')
file12=open(pwd3,'w')
list1=[]
for gseq in SeqIO.parse(pwd11,'fasta'):
	cp=gseq.description.split('\t')[1].split('-')[1]
	if len(gseq.seq)>400 and float(cp)>=5:
		list1.append(gseq)

for gseq in SeqIO.parse(pwd1,'fasta'):
	i1=gseq
l='>'+i1.description+"\n"+i1.seq+"\n"
#file1.writelines(l)
pwda='temp1.fa'
ft=open(pwda,'w')
ft.writelines(l)
ft.close()
j0=0
list11=[]
for i2 in list1:
	id3=i2.id.split('.unit-')[0]
	if True:#id3==i1.id:
	#	l='>'+i2.description+"\t+\n"+i2.seq+"\n"
	#	file1.writelines(l)
	#else:
	#	cp=i2.description.split('\t')[1].split('-')[1]
	#	if float(cp)<5:
	#		continue
	#	else:
	#		print (cp)
		print (i2.id,j0)
		j0+=1
		pwdb="temp2.fa"
		l='>'+i2.description+"\n"+str(i2.seq).replace("-", '')+"\n"
		ft=open(pwdb,'w')
		ft.writelines(l)
		ft.close()

		commandline='python ./mummer_2.py --seq2 '+pwda+' --seq1 '+pwdb+' --name temp3'
		os.system(commandline)

		file2=open('temp3.delta.filter.coords','r')
		list2=[]
		flag=1
		for row in file2:
			row1=row.rstrip().split('\t')
			if float(row1[6])>=55.0:
				print (row,int(row1[2])==1,int(row1[2])==int(row1[8]),int(row1[3])==1,int(row1[3])==int(row1[8]))
				if int(row1[2])==1 or int(row1[2])==int(row1[8]) or int(row1[3])==1 or int(row1[3])==int(row1[8]):
					if list2==[]:
						list2=row1
					else:
						if int(list2[5])<=int(row1[5]):
							list2=row1
					flag=1
		if list2==[]:
			file2.seek(0,0)
			for row in file2:
				row1=row.rstrip().split('\t')
				if float(row1[6])>=55.0:
					if int(row1[0])==1 or int(row1[0])==int(row1[7]) or int(row1[1])==1 or int(row1[1])==int(row1[7]):
						if list2==[]:
							list2=row1
						else:
							if int(list2[5])<=int(row1[5]):
								list2=row1
						flag=0
		if list2==[] :#or flag==0:
			print (list2)
			print (i2.seq)
			print (i1.seq)
			file2.seek(0,0)
			for row in file2:
				row1=row.rstrip().split('\t')
				if list2==[]:
					list2=row1
				else:
					if int(list2[5])<=int(row1[5]):
						list2=row1
			print (list2)
	#		sys.exit()
		#if 'GWHBWDJ00000013.1' in row:
		#	print (list2)
		#	sys.exit()
		file2.close()
		commandline='rm temp2*'
		os.system(commandline)
		commandline='rm temp3*'
		os.system(commandline)
		if list2!=[]:
			print (list2)
			s=int(list2[2])-1
			e=int(list2[3])
			seq1=i2.seq
			s1=int(list2[0])-1
			e1=int(list2[1])
			d="+"
			if s>e:
				seq1=i2.seq.reverse_complement()
				s=int(list2[3])-1
				e=int(list2[2])
				e1=int(list2[7])-int(list2[0])+1
				s1=int(list2[7])-int(list2[1])
				d="-"
			d1=s-0
			d2=int(list2[8])-e
			dt=min(s-0,int(list2[8])-e)
			print (d1,d2,d1<=d2,dt)
			if d1<=d2: #s=0
				#seq2=seq1[s1:]+seq1[:s1]
				seq2=seq1[s1-dt:]+seq1[:s1-dt]
				if len(seq2)!=len(seq1):
					print ('error')
					sys.exit()
			else:
				print (e1-dt)
				#seq2=seq1[e1:]+seq1[:e1]
				seq2=seq1[e1+dt:]+seq1[:e1+dt]
				if len(seq2)!=len(seq1):
					print ('error')
					sys.exit()
			l='>'+i2.description+"\t"+d+"\n"+seq2+"\n"
			ft=open(pwdb,'w')
			ft.writelines(l)
			ft.close()
			commandline='python ./mummer_1.py --seq2 '+pwda+' --seq1 '+pwdb+' --name temp3'
			os.system(commandline)
			file2=open('temp3.delta.filter.coords','r')
			list2=[]
			tag=0
			print('after')
			
			for row in file2:
				print (row)
				#list2.append(row)
				row1=row.rstrip().split('\t')
				if float(row1[9])>=95.0 and float(row1[10])>=95.0:
					tag=1
					list2.append(row)
			commandline='rm temp2*'
			os.system(commandline)
			commandline='rm temp3*'
			os.system(commandline)
			print (i2.seq)
			print (seq2,'F')
			print (i1.seq,'F')
			#if flag==0:
			#	sys.exit()
			if len(list2)!=1 and tag==0:
				print (list2)
				file12.writelines(l)
			#	file1.writelines(l)
			elif tag==1:
				seq3=str(seq2).replace("-", '')
				l='>'+i2.description+row.rstrip()+"\t"+d+"\n"+seq3+"\n"
				file1.writelines(l)
			#if 'GWHBWDJ00000018.1' in row:
			#a	sys.exit()
		#	if d1>d2:
		#		sys.exit()	
		else:
			list11.append(i2.seq)
#print (list11)
#for i in list11:
#	print (i)
file1.close()
file12.close()
print ('welldone-all')
