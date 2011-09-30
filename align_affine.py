"""This is a Python script for aligning two aminoacid sequences 
using the Needleman-Wunch algorithm for global alignment, and the
Smith-Waterman algorithm for local alignment. It uses an affine 
gap model.

Not yet finished. Still having some problems with the Iy and Ix
matrices. And no backtracking yet. And no local alignment. Yet. 

Author: Paula Petcu
Created: 28 September 2011
Licensed under: GNU General Public License v2
"""

import os
import sys
from pprint import pprint
import pickle
from sets import Set

MININF = float("-inf") # minus infinity in Pyhton

scoringMatrices = ["blosum50"] #Currently only supporting BLOSUM50
aminoAcids = "ARNDCQEGHILKMFPSTWYV"
trace = {-1:'-',0:'\\',1:'|',2:'-',3:' '}

def readScoringMatrix(matrixName):
	"""Read the scoring matrix (either from a pickle file, if it exists, or from 
	a text file - if none exists, returns None). 
	"""
	
	matrix = None
	if os.path.isfile(matrixName+".p"):
		pckl_file = open(matrixName+".p","rb")
		matrix = pickle.load(pckl_file)
		pckl_file.close()
	elif os.path.isfile(matrixName+".txt"):
		lines = open(matrixName+".txt","r").readlines()
		matrix = []
		for line in lines:
			matrix.append(line.rstrip().split("\t"))
		pckl_file = open(matrixName+".p","wb")
		pickle.dump(matrix,pckl_file)
		pckl_file.close()
	return matrix
	
def areArgumentsValid(seqA,seqB,gap,expand,match,mismatch,matrixName):
	"""Verify the validity of the arguments
	"""
	
	if( (match!=None and mismatch!=None and matrixName!=None) or\
	 (match==None and mismatch==None and matrixName==None)):
		print "Too many or no scoring values given!"
		return False
	if matrixName and matrixName not in scoringMatrices:
		print "Do not have the values for this matrix."
		return False
	try:
		int(gap)
		int(expand)
	except:
		print "The gap and expand values should be positive integers."
		return False
	if match!=None and mismatch!=None:
		try:
			int(match)
			int(mismatch)
		except:
			print "The match and mismatch values should be integers."
			return False		
	if not set(aminoAcids).issuperset(set(seqA)) or not set(aminoAcids).issuperset(set(seqB)):
		print "The sequences must be composed of aminoacids ("+aminoAcids+")."
		return False
	return True
	
def needle(seqA,seqB,gap,extend,match=None,mismatch=None,matrixName=None):
	"""Apply Needleman-Wunch on two AA sequences (global alignment)
	Using the affine gap model.
	"""
	
	#Verify the arguments
	if not areArgumentsValid(seqA,seqB,gap,extend,match,mismatch,matrixName):
		return
	gap=int(gap)
	extend=int(extend)
	if match!=None and mismatch!=None:
		match = int(match)
		mismatch = int(mismatch)
	
	print "Aligning sequence "+seqA+" with sequence "+seqB+" using Needleman-Wunch:"
	
	#Compute the matrix of pair scores
	if matrixName==None:
		pairScores = createMatrixOfPairScores(seqA,seqB,match,mismatch)
	else:
		blosum50 = readScoringMatrix(matrixName)
		pairScores = createMatrixOfPairScores(seqA,seqB,matrix=blosum50)
	print "Here is the matrix of pair scores for the two sequences: "
	pprint(pairScores)
	
	#Generate the dynamic programming matrix (global align., affine gap model)
	Iy = []
	M = []
	Ix = []
	# - initialize M, Ix, and Iy for first row
	Iy.append([MININF]*(len(seqA)+1))
	M.append([0])
	M[0] += [MININF]*len(seqA)
	Ix.append([MININF])
	Ix[0] += range(-gap,-gap-extend*len(seqA),-extend)
	# - initialize tracking matrices
	trackIy = []
	trackM = []
	trackIx = []
	trackIy.append([trace[3]]*(len(seqA)+1)) # put ' ' on first Iy row
	trackM.append([trace[3]]*(len(seqA)+1)) # put ' ' on first M row
	trackIx.append([' ']+['--']*(len(seqA))) # put '--' on first Ix row
	for i,aB in enumerate(seqB):
		Iy.append([])
		M.append([])
		Ix.append([])
		# - initialize M, Ix, and Iy for first column
		Iy[i+1].append(-gap-extend*i)
		M[i+1].append(MININF)
		Ix[i+1].append(MININF)
		# - initialize tracking matrices for first column
		trackIy.append([])
		trackIy[i+1].append('||') # put '||' on first Iy column
		trackM.append([])
		trackM[i+1].append(trace[3]) # put ' ' on first M column
		trackIx.append([])
		trackIx[i+1].append(trace[3]) # put ' ' on first Ix column
		for j,aA in enumerate(seqA):
			maxArgs_Ix = [M[i][j+1]-gap,\
						Ix[i][j+1]-extend]
			maxVal_Ix = max(maxArgs_Ix)
			maxArgs_Iy = [M[i+1][j]-gap,\
						Iy[i+1][j]-extend]
			maxVal_Iy = max(maxArgs_Iy)
			maxArgs = [M[i][j] + pairScores[i][j],\
				Ix[i][j] + pairScores[i][j],\
				Iy[i][j] + pairScores[i][j]]
			maxVal = max(maxArgs)
			Iy[i+1].append(maxVal_Iy)
			M[i+1].append(maxVal)
			Ix[i+1].append(maxVal_Ix)
			# - update tracking matrices
			if maxArgs_Iy.index(maxVal_Iy) == 0: trackIy[i+1].append('/')
			else: trackIy[i+1].append('-')
			if maxArgs.index(maxVal) == 0: trackM[i+1].append('\\')
			elif maxArgs.index(maxVal) == 1: trackM[i+1].append('-')
			else: trackM[i+1].append('|')
			if maxArgs_Ix.index(maxVal_Ix) == 0: trackIx[i+1].append('|')
			else: trackIx[i+1].append('\\')
	
	print "Here is the global dynamic programming matrix for the two sequences: "
	print "Iy="
	pprint(Iy)
	print "M="
	pprint(M)
	print "Ix="
	pprint(Ix)
	
	print "Here is the global dynamic programming traceback for the two sequences: "
	print "trackIy="
	pprint(trackIy)
	print "trackM="
	pprint(trackM)
	print "trackIx="
	pprint(trackIx)
	
	#Backtracking
	# Start from bottom right corner
	loc = (len(seqB),len(seqA))
	maxArgs = [Iy[loc[0]][loc[1]], M[loc[0]][loc[1]], Ix[loc[0]][loc[1]]]
	score = max(maxArgs)
	which = maxArgs.index(score) 
	alignment =[[],[]]
	# Go towards the top left corner
	while loc != (0,0):
		#print loc, which
		#print alignment
		if which==0:
			direction = trackIy[loc[0]][loc[1]]
			if direction=='-' or direction=='/': 
				if direction=='-': which = 0 #from Iy
				elif direction=='/': which = 1 #from M
				alignment[0].insert(0,seqA[loc[1]-1])
				alignment[1].insert(0,'-')
				loc = (loc[0],loc[1]-1)
			elif direction=='||': #margin hit, go up
				which = 1
				alignment[0].insert(0,'-')
				alignment[1].insert(0,seqB[loc[0]-1])
				loc = (loc[0]-1,loc[1])
		elif which==2:
			direction = trackIx[loc[0]][loc[1]]
			if direction=='|' or direction=='\\':
				if direction=='|': which = 1 #from M
				elif direction=='\\': which = 2 #from Ix
				alignment[0].insert(0,'-')
				alignment[1].insert(0,seqB[loc[0]-1])
				loc = (loc[0]-1,loc[1])
			elif direction=='--': # margin hit, go left
				which = 2
				alignment[0].insert(0,seqA[loc[1]-1])
				alignment[1].insert(0,'-')
				loc = (loc[0],loc[1]-1)
		elif which==1:
			direction = trackM[loc[0]][loc[1]]
			alignment[0].insert(0,seqA[loc[1]-1])
			alignment[1].insert(0,seqB[loc[0]-1])
			loc = (loc[0]-1,loc[1]-1)
			if direction=='|': which = 0 #from Iy
			elif direction=='\\': which = 1 #from M
			elif direction=='-': which = 2 #from Ix
	print "Here is a global alignment:"
	#pprint(alignment)
	for i in range(len(alignment)): print ''.join(alignment[i])
	
	return

def createMatrixOfPairScores(seqA,seqB,match=None,mismatch=None,matrix=None):
	"""Create the matrix of pair scores for two sequences, based either on match/mismatch
	values or on a scoring matrix.
	"""
	
	pairScores = []
	if matrix==None:
		for i,aB in enumerate(seqB):
			pairScores.append([])
			for j,aA in enumerate(seqA):
				if aA == aB:
					pairScores[i].append(match)
				else:
					pairScores[i].append(mismatch)
	else:
		for i,aB in enumerate(seqB):
			pairScores.append([])
			for j,aA in enumerate(seqA):
				pairScores[i].append(int(matrix[matrix[0].index(aA)][matrix[0].index(aB)]))			
	return pairScores

if len(sys.argv)!=1 and len(sys.argv)!=7 and len(sys.argv)!=8:
	print "Usage: python align_affine.py"
	print "Usage: python align.py global blossum50 gap extend HEAGAWGHEE PAWHEAE"
	print "Usage: python align.py global match mismatch gap extend HEAGAWGHEE PAWHEAE"
	exit(0)
elif len(sys.argv)==1:
	#Default behaviour
	# first example
	seqA = "HEAGAWGHEE"
	seqB = "PAWHEAE"
	#needle(seqA,seqB,gap=2,match=1,mismatch=-2)
	needle(seqA,seqB,gap=10,extend=1,matrixName="blosum50")
	# a second example
	seqA = "GHGKKVADALTN"
	seqB = "GHKRLLT"
	needle(seqA,seqB,gap=10,extend=1,matrixName="blosum50")
	# a third example
	seqA = "VLSPADK"
	seqB = "HLAESK"
	needle(seqA,seqB,gap=12,extend=2,matrixName="blosum50")
'''elif len(sys.argv)==6 or len(sys.argv)==7:
	#Alignment according to the received arguments
	if sys.argv[1]=="global":
		if len(sys.argv)==6:
			needle(sys.argv[4],sys.argv[5],gap=sys.argv[3],matrixName=sys.argv[2])
		else:
			needle(sys.argv[5],sys.argv[6],gap=sys.argv[4],match=sys.argv[2],mismatch=sys.argv[3])
	elif sys.argv[1]=="local":
		if len(sys.argv)==6:
			water(sys.argv[4],sys.argv[5],gap=sys.argv[3],matrixName=sys.argv[2])
		else:
			water(sys.argv[5],sys.argv[6],gap=sys.argv[4],match=sys.argv[2],mismatch=sys.argv[3])
	else:
		print "Usage: python align.py global blossum50 gap HEAGAWGHEE PAWHEAE"
		print sys.argv[1]+" is not a valid arg. Use global or local."
		exit(0)
'''
