import re
pwd1='Jack.uniq.grouped.txt'
file1=open(pwd1,'r')
dict1={}
for row in file1:
	row1=row.rstrip().split('\t')
	ID=row1[-1].split(';')
	dict1[row1[1]]=ID
file1.close()
pwd1='Jack.Type.grouped.txt'
file1=open(pwd1,'r')
pwd2='Jack.Type.grouped.all.txt'
file2=open(pwd2,'w')
for row in file1:
	row1=row.rstrip().split('\t')
	ID=row1[-1].split(';')
	listD1=[]
	for i in ID:
		Dt=dict1[i]
		listD1=listD1+Dt
	l="\t".join(row1[:-1])+"\t"+';'.join(listD1)+"\n"
	file2.writelines(l)
file1.close()
file2.close()
print ('welldone')
