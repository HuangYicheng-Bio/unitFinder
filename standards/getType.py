import getopt
import sys
import os
import Bio
from Bio import SeqIO
from mummer_2 import mummer2
from completeGm1 import completeGm1
from completeGm4 import completeGm4
from getHORBlock import getHORBlock
from getDistenseBlock import getDistenseBlock
from aln_1 import aln1
from aln_4 import aln4
from uniqStandards import uniqStandards
def usage():
	print ("--in *.Consensus.fa")
	print ("--unit *.Units.fa")
	print ("--type1 centGm1.fa")
	print ("--type4 centGm4.fa")
	print ("--name ID")
	print ("-h help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "h",["in=","name=","unit=","type1=","type4=","RM=","con=","cov="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()

	seq1=''
	units=''
	type1=''
	type4=''
	name=''
	RM=''
	con=''
	cov=''
	if len(opts)==0:
		print ("WHERE IS PARAMETER!!!")
		print ("Use '--h' for same information")
		sys.exit()
	for opt,value in opts:
		if opt in ("--in"):
			seq1=value
		if opt in ("--unit"):
			units=value
		if opt in ("--type1"):
			type1=value
		if opt in ("--type4"):
			type4=value
		elif opt in ("--name"):
			name=value
		elif opt in ("--RM"):
			RM=value
		elif opt in ("--con"):
			con=value
		elif opt in ("--cov"):
			cov=value
		elif opt in ("-h"):
			usage()
			sys.exit()
	#getTypeFromMerge(seq1,seq2,type1,name)
	return seq1,units,type1,type4,name,RM,con,cov
seq1,units,type1,type4,name,RM,con,cov=getoptions()
commond="mkdir Aln"
os.system(commond)
if os.path.exists('./Aln/'+name+'-CentGm-1.delta.filter.coords')!=True or os.path.getsize('./Aln/'+name+'-CentGm-1.delta.filter.coords')==0:
	mummer2(type1,units,"./Aln/"+name+"-CentGm-1")
if os.path.exists('./Aln/'+name+'-CentGm-4.delta.filter.coords')!=True or os.path.getsize('./Aln/'+name+'-CentGm-4.delta.filter.coords')==0:
	mummer2(type4,units,"./Aln/"+name+"-CentGm-4")

commond="mkdir Standards"
os.system(commond)
if os.path.exists('./Standards/CentGm-1-'+name+'.complete.fa')!=True or os.path.getsize('./Standards/CentGm-1-'+name+'.complete.fa')==0:
	#sys.exit()
	completeGm1(name,units,type1,seq1,'./Aln/'+name+'-CentGm-1.delta.filter.coords')
	uniqStandards('./Standards/CentGm-1-'+name+'.complete.fa','./Standards/CentGm-1-'+name+'.complete.uniq.fa')
#sys.exit()
if os.path.exists('./Standards/CentGm-4-'+name+'.complete.fa')!=True or os.path.getsize('./Standards/CentGm-4-'+name+'.complete.fa')==0:
	#sys.exit()
	completeGm4(name,units,type4,seq1,'./Aln/'+name+'-CentGm-4.delta.filter.coords')
	uniqStandards('./Standards/CentGm-4-'+name+'.complete.fa','./Standards/CentGm-4-'+name+'.complete.uniq.fa')

commond="cat ./Standards/CentGm-*.complete.uniq.fa >CentGm-units.fa"
os.system(commond)
for gseq in SeqIO.parse(seq1,'fasta'):
	loc=gseq.id.split('-')
	start=int(loc[-2])
	end=int(loc[-1])
	ID1="-".join(loc[:-2])
#commond="cat "+RM+" | grep '"+ID1+"' | awk '{if($7>="+str(start)+" && $6<="+str(end)+"){print}}' | awk '{if ($12==\"(0)\"|| $14==\"(0)\"){print}}' >"+name+".RM.1.txt"
commond="cat "+RM+" | grep '"+ID1+"' | awk '{if($7>="+str(start)+" && $6<="+str(end)+"){print}}'  >"+name+".RM.1.txt"
print (commond)
os.system(commond)

#if os.path.exists(name+'.block.units.fa')!=True or os.path.getsize(name+'.block.units.fa')==0:
	#getHORBlock(con,"CentGm-units.fa",name+".RM.1.txt",name,start,end)
	#getHORBlock(con,units,name+".RM.1.txt",name,start,end)
getHORBlock(con,'./Standards/CentGm-1-'+name+'.complete.uniq.fa','./Standards/CentGm-4-'+name+'.complete.uniq.fa',name+".RM.1.txt",name,start,end)
getDistenseBlock(name+'.block.units.fa',name+".RM.1.txt",name+".block",start,end,ID1,"CentGm-units.fa",cov,'./Standards/CentGm-1-'+name+'.complete.uniq.fa','./Standards/CentGm-4-'+name+'.complete.uniq.fa')
#getDistenseBlock(name+'.block.units.fa',name+".RM.1.txt",name+".block",start,end,ID1,units,cov)
if os.path.exists('CentGm-1-'+name+'.type.fa')!=True or os.path.getsize('CentGm-1-'+name+'.type.fa')==0:
	aln1('./Standards/CentGm-1-'+name+'.complete.fa','CentGm-1-'+name+'.type.fa',name)
if os.path.exists('CentGm-4-'+name+'.type.fa')!=True or os.path.getsize('CentGm-4-'+name+'.type.fa')==0:
	aln4('./Standards/CentGm-4-'+name+'.complete.fa','CentGm-4-'+name+'.type.fa',name)
