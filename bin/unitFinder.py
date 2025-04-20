import getopt
import sys
import os
from mummer_cen import mummer
from mergeAln import mergeAln
from plotRepeatGroup import plotRepeatGroup
from groupTRF import groupTRF
from getseq import getseq
from secondary import secondary
from length import length
from uniqUnits import uniqUnits
from plotRepeatBlock import plotRepeatBlock
from typeMerge import typeMerge
from getUniqBlock import getUniqBlock
from getTypeFromMerge import getTypeFromMerge
def usage():
	print ("--seq *.fa")
	print ("--name ID")
	print ('--cen *.fa(optional)')
	print ("--flag 0|1 (default 0)\n\t\t0:find potential centromere region.\n\t\t1:find repeats with whole sequence.")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["name=","seq=","cen=","flag="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	name=''
	cen=''
	flag=0
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--seq"):
			seq1=value
		elif opt in ("--name"):
			name=value
		elif opt in ("--cen"):
			cen=value
		elif opt in ("--flag"):
			flag=int(value)
		elif opt in ("-h"):
			usage()
			sys.exit()
		
	return name,seq1,cen,flag

def unitFinder(start,rowoutpwd,name):
	commandline='mkdir findTandemRepeats'
	if os.path.exists('findTandemRepeats')!=True:
		os.system(commandline)
	
	if os.path.exists(name+'.Units.fa')!=True or os.path.getsize(name+'.Units.fa')==0:
		secondary(rowoutpwd,name)
		commandline='mv round* ./findTandemRepeats/'
		os.system(commandline)
		commandline='mv unit_round* ./findTandemRepeats/'
		os.system(commandline)
		commandline='mv '+name+'.loci-* ./findTandemRepeats/'
		os.system(commandline)
	
	length(name+'.Units.fa',name+'.Units.length')
	plotRepeatBlock(name+'.Units.fa',name+'.Consensus.fa',name+'.units',start)
	uniqUnits(name+'.Units.fa',name+'.Consensus.fa',name)
	length(name+'.Units.uniq.fa',name+'.Units.uniq.length')

	plotRepeatBlock(name+'.Units.uniq.fa',name+'.Consensus.uniq.fa',name+'.uniq.units',start)
	#getUniqBlock(name+'.Consensus.fa',name+'.Units.uniq.fa',name+'.units.stat.txt',name)
	#plotRepeatBlock(name+'.block.units.fa',name+'.block.fa',name+'.block.units',start)
	#length(name+'.block.units.fa',name+'.block.units.length')
	Grouppwd=name+'.Consensus.uniq.grouped.txt'
	if os.path.exists(Grouppwd)!=True or os.path.getsize(Grouppwd)==0:
		#print (os.path.exists(Grouppwd)!=True, os.path.getsize(Grouppwd)==0)
		#sys.exit()
		typeMerge(name+'.Consensus.uniq.fa',name+'.Consensus.uniq.grouped.txt',name)
	getTypeFromMerge(name+'.Consensus.uniq.fa', name+'.Units.uniq.fa',name+'.Consensus.uniq.grouped.txt',name)
	plotRepeatGroup(name+'.Type.units.fa',name+'.Consensus.uniq.grouped.txt',start,name+'.Type.units')


name,seq1,cen,flag=getoptions()
if flag==0:
	if cen!='':
		Aln=name+".delta.filter.coords"
		if os.path.exists(Aln)!=True or os.path.getsize(Aln)==0:
			Aln=mummer(cen,seq1,name)
		start,end=mergeAln(Aln)
		print ('Location:',start,end)
		rowoutpwd=name+".loci.fa"
		getseq(seq1,int(start),int(end),rowoutpwd)
	else:
		commandline='trf '+seq1+' 2 7 7 80 10 50 500 -f -d -m'
		nametrf=seq1.split('/')[-1]
		pwdtrf=nametrf+'.2.7.7.80.10.50.500.dat'
		if os.path.exists(pwdtrf)!=True or os.path.getsize(pwdtrf)==0:
			pwdtrf1='./coarsePosition/'+pwdtrf
			if os.path.exists(pwdtrf1)!=True or os.path.getsize(pwdtrf1)==0:
				os.system(commandline)
			else:
				pwdtrf=pwdtrf1
		start,end=groupTRF(name,pwdtrf)
		print ('Location:',start,end)
		rowoutpwd=name+".loci.fa"
		getseq(seq1,int(start),int(end),rowoutpwd)
		if 'coarsePosition' not in pwdtrf:
			commandline='mkdir coarsePosition'
			os.system(commandline)
			commandline='mv *2.7.7.80.10.50.500* ./coarsePosition'
			os.system(commandline)
else:
	start=0
	rowoutpwd=seq1
unitFinder(start,rowoutpwd,name)
