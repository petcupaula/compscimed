"""This is a Python script for reading the data file that can be downloaded 
from 23andMe. It reads the file and transforms it into a pickled file.

Author: Paula Petcu
Created: 24 September 2011
Licensed under: GNU General Public License v2
"""

import cPickle
import sys

def readTxtFile(genomeFileName,snpFileName):
	"""Read a 23andMe txt file with TAB-separated snp data into a list and 
	dump the structure using pickle.
	"""
	
	genomeFile = open(genomeFileName,"r")
	allLines = genomeFile.readlines()
	genomeFile.close()

	snps = [["rsid","chromosome","position","genotype"]]
	for i,line in enumerate(allLines):
		#if i%1000==0: print i
		if not line[0]=="#":
			info = line.rstrip().split('\t')
			snps.append(info)
	
	snpFile = open(snpFileName,"wb")
	cPickle.dump(snps,snpFile)
	snpFile.close()
	
	print "Output generated in file "+snpFileName
	return snpFileName
	
if len(sys.argv)!=2 and len(sys.argv)!=3:
	print "Usage: txt2p.py genome_pa.txt"
	print "Usage: txt2p.py genome_pa.txt genome_pa.p"
	exit(0)
elif len(sys.argv)==2:
	readTxtFile(genomeFileName=sys.argv[1],snpFileName=sys.argv[1].replace(".txt",".p"))
else:
	readTxtFile(genomeFileName=sys.argv[1],snpFileName=sys.argv[2])
	