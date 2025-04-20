import re
import statistics
import getopt
import sys
import numpy as np
def usage():
	print ("--trf *.dat")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["name=","trf="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq2=''
	name=''
	trf=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--trf"):
			trf=value
		elif opt in ("--name"):
			name=value

		elif opt in ("-h"):
			usage()
			sys.exit()
	return name,trf

def groupTRF(name,trf):

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
		for j in lt:
			j.append(i)
		lt1=lt1+lt
	lTRF=[]

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
					break
				ltn+=1
			if flag==0:
				lTRF.append([a1,b1,c1])
	print ('cluster number:',len(lTRF))

	l=[]
	j0=1
	lmax=[]
	for i in lTRF:
		a,b,c=i
		l.append([c,a,b,[str(j0)]])
		if lmax==[]:
			lmax=[[c,a,b,[str(j0)]]]
		else:
			if c>lmax[0][0]:
				lmax=[[c,a,b,[str(j0)]]]
		j0+=1
	
	start=lmax[0][1]
	end=lmax[0][2]
	
	print ('welldone')
	return start,end
