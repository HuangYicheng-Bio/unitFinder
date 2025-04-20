import os
import sys
import getopt

def usage():
	print ("--seq1 seq1")
	print ("--seq2 seq2")
	print ("--name output")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["seq1=","seq2=","name="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	seq2=''
	name=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--seq1"):
			seq1=value
		elif opt in ("--seq2"):
			seq2=value
		
		elif opt in ("--name"):
			name=value
		
		elif opt in ("-h"):
			usage()
			sys.exit()
	return seq1,seq2,name

def mummer(seq1,seq2,name):
	mr=name
	commendline="nucmer -c 5 -l 5 -p "+mr+" "+seq1+" "+seq2
	#print (commendline)
	os.system(commendline)
	commendline="delta-filter -m "+mr+".delta > "+mr+".delta.filter"
	#print (commendline)
	os.system(commendline)
	commendline="show-coords -TrHcl "+mr+".delta.filter > "+mr+".delta.filter.coords"
	#print (commendline)
	os.system(commendline)
	mr=mr+".delta.filter"
	return mr

def show(pwd):
	pwd=pwd+".coords"
	ft=open(pwd,'r')
	l=ft.readline()
	while l:
		print (l)
		l=ft.readline()
	ft.close()
	return 0

def mummerplot(seq1,seq2,name,mr):
	cl="mummerplot -t png --prefix="+name+" "+mr+" -R "+seq1+" -Q "+seq2+" --filter --layout"
	print (cl)
	os.system(cl)
	return 0

seq1,seq2,name=getoptions()
r=mummer(seq1,seq2,name)
#show(r)
#mummerplot(seq1,seq2,name,r)
#print ("well done")

