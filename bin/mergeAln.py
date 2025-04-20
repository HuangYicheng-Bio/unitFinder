import re
import sys

def merge(list1,mercy):
	list2=sorted(list1,key=lambda x:x[0])
	j=0
	while j<len(list2):
		j1=j+1
		while j1<len(list2):
			a,b,n1=list2[j]
			c,d,n2=list2[j1]
			t1=-1
			t2=-1
			if c>=a and c-mercy<b:
				t1=a
				if b>d:
					t2=b
				else:
					t2=d
			else:
				if a>=c and a-mercy<d:
					t1=c
					if d>b:
						t2=d
					else:
						t2=b
				
			if t1!=-1 and t2!=-1:
				if t1<0 or t2<0:
					sys.exit()
				list2[j]=[t1,t2,n1+n2]
				list2.pop(j1)
				j1-=1
			j1+=1
		j+=1
	list3=[]
	listLen=0
	for i in list2:
		s,e,n=i
		if s!=e:
			list3.append(i)
			listLen=listLen+(e-s+1)
	return list3,listLen

def mergeAln(Aln):
	file1=open(Aln,'r')
	list1=[]
	for row in file1:
		row1=row.rstrip().split('\t')
		a=int(row1[2])
		b=int(row1[3])
		c=row1[-2].split('_Repeat_')[0]
		if a>b:
			a=int(row1[3])
			b=int(row1[2])
		list1.append([a,b,[c]])
	list2,len1=merge(list1,10000)
	list22=[]
	for i in list2:
		a=i[1]-i[0]
		#print (i[0],i[1],len(set(i[2])),a)
		if len(set(i[2]))>=10 or a>=100000:#and a>=100000:
			print (i[0],i[1],len(set(i[2])),a)
			list22.append(i)
	file1.close()
	list3,len1=merge(list22,100000)
	max1=-1
	lams=[]
	lame=[]
	lamm=[]
	for i in list3:
		a=i[1]-i[0]
		if len(set(i[2]))>=10:
			print (i[0],i[1],len(set(i[2])),a)
			lams.append(i[0])
			lame.append(i[1])
		if a>max1:
			max1=a
			lamm=i
	if lams!=[]:
		return min(lams),max(lame)
	else:
		return lamm[0],lamm[1]
