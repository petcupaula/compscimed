"""This is a Python script that reads the pickled data file downloaded from
23andMe, and extracts chromosome 15 data.

Run txt2p.py first.

Author: Paula Petcu
Created: 1 October 2011
Licensed under: GNU General Public License v2
"""

import cPickle
import sys

def extractChr15(snpFileName,chr15FileName):
	""" Save the data for chromosome 15 in a file. snpFileName is the
	file with the entire list of personal snp data
	"""
	
	# Read SNP pickle file (first row is ["rsid","chromosome","position","genotype"])
	snpFile = open(snpFileName,"rb")
	snps = cPickle.load(snpFile)
	snpFile.close()
	
	# Select data for chromosome 15
	chr15 = [snps[0]] # first line are the field names
	for line in snps[1:]:
		if line[1] == "15":
			chr15.append(line)
	
	# Dump chromosome 15 data
	chr15File = open(chr15FileName,"wb")
	cPickle.dump(chr15,chr15File)
	chr15File.close()
	
	print "Output generated in file "+chr15FileName
	return chr15FileName

if len(sys.argv)!=2 and len(sys.argv)!=3:
	print "Usage: extract15.py genome_pa_snp.p"
	print "Usage: extract15.py genome_pa_snp.p chr15.p"
	exit(0)
elif len(sys.argv)==2:
	extractChr15(snpFileName=sys.argv[1],chr15FileName="chr15.p")
else:
	extractChr15(snpFileName=sys.argv[1],chr15FileName=sys.argv[2])
