"""This is a Python script that prints some statistics about the data
contained in a pickled 23andMe data file (i.e. total number of SNPs,
number of SNPs per chromosome, etc.) 

Run txt2p.py first.

Author: Paula Petcu
Created: 24 September 2011
Licensed under: GNU General Public License v2
"""

import cPickle
import os
import sys

def printStatistics(snpFileName):
	"""Read the dumped SNP data and print some statistics
	"""
	
	snpFile = open(snpFileName,"rb")
	snps = cPickle.load(snpFile)
	snpFile.close()
	
	print "# Number of SNPs: "
	print str(len(snps)-1) # first row is column title
	
	RS = 0
	I = 0
	NC = 0
	snpsPerChrom = {'X':[0,0,0,0],'Y':[0,0,0,0],'MT':[0,0,0,0]}
	for i in range(1,23): snpsPerChrom[str(i)] = [0,0,0,0]
	for line in snps[1:]:
		snpsPerChrom[line[1]][0] += 1
		if line[0].startswith("rs"): 
			RS += 1
			snpsPerChrom[line[1]][1] += 1
		elif line[0].startswith("i"): 
			I += 1
			snpsPerChrom[line[1]][2] += 1
		if line[3] == '--':
			NC += 1
			snpsPerChrom[line[1]][3] += 1
	
	print "# Number of rs identifiers: "
	print str(RS)
	print "# Number of internal identifiers: "
	print str(I)
	print "# Number of no-calls: "
	print str(NC)
	print "# Number of SNPs per chromosome (total number of snps, number of rsid, \
	number of iid, number of no-call genotypes): "
	#print str(snpsPerChrom)
	for elem in snpsPerChrom:
		print elem+" "+" ".join(str(snpsPerChrom[elem])[1:-1].split(","))

if len(sys.argv)!=2:
	print "Usage: python snpstats.py snps.p"
	exit(0)
else:
	if not sys.argv[1].endswith('.p'):
		print "Not a valid filetype. Run 'python txt2p.py genome_you.txt' first."
		exit(0)
	snpFileName = sys.argv[1]
	printStatistics(snpFileName)
