"""This is a Python script for aligning two aminoacid sequences 
using the Needleman-Wunch algorithm for global alignment, and the
Smith-Waterman algorithm for local alignment. It uses a linear 
gap model.

Author: Paula Petcu
Created: 24 September 2011
Licensed under: GNU General Public License v2
"""

import os
import sys
from pprint import pprint
import pickle
from sets import Set

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
	
def areArgumentsValid(seqA,seqB,gap,match,mismatch,matrixName):
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
	except:
		print "The gap value should be a positive integer."
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
	
def needle(seqA,seqB,gap,match=None,mismatch=None,matrixName=None):
	"""Apply Needleman-Wunch on two AA sequences (global alignment)
	Using the linear gap model.
	"""
	
	#Verify the arguments
	if not areArgumentsValid(seqA,seqB,gap,match,mismatch,matrixName):
		return
	gap=int(gap)
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
	
	#Generate the dynamic programming matrix (global align., linear gap model)
	# - initialize the first row of the dp matrix using the gap value 
	dpMatrix = []
	dpMatrix.append(range(0,-gap*(len(seqA)+1),-gap))
	dpMatrixTrack = []
	dpMatrixTrack.append([trace[2]]*(len(seqA)+1)) # - put '-' on first row
	for i,aB in enumerate(seqB):
		dpMatrix.append([])
		dpMatrix[i+1].append(-gap*(i+1))
		dpMatrixTrack.append([])
		dpMatrixTrack[i+1].append(trace[1]) # - put '|' on first column
		for j,aA in enumerate(seqA):
			maxArgs = [dpMatrix[i][j] + pairScores[i][j],\
				dpMatrix[i][j+1] - gap,\
				dpMatrix[i+1][j] - gap]
			maxVal = max(maxArgs)
			dpMatrix[i+1].append(maxVal)
			dpMatrixTrack[i+1].append(trace[maxArgs.index(maxVal)])
	print "Here is the global dynamic programming matrix for the two sequences: "
	pprint(dpMatrix)
	print "Here is the global dynamic programming traceback for the two sequences: "
	pprint(dpMatrixTrack)
	
	#Backtracking
	# Start from bottom right corner
	loc = (len(seqB),len(seqA))
	alignment =[[],[]]
	# Go towards the top left corner
	while loc != (0,0):
		#print alignment
		direction = dpMatrixTrack[loc[0]][loc[1]]
		if direction=='-':
			alignment[0].insert(0,seqA[loc[1]-1])
			alignment[1].insert(0,'-')
			loc = (loc[0],loc[1]-1)
		elif direction=='|':
			alignment[0].insert(0,'-')
			alignment[1].insert(0,seqB[loc[0]-1])
			loc = (loc[0]-1,loc[1])
		elif direction=='\\':
			alignment[0].insert(0,seqA[loc[1]-1])
			alignment[1].insert(0,seqB[loc[0]-1])
			loc = (loc[0]-1,loc[1]-1)		
	print "Here is a global alignment:"
	#pprint(alignment)
	for i in range(len(alignment)): print ''.join(alignment[i])
	
	return

def water(seqA,seqB,gap,match=None,mismatch=None,matrixName=None):
	"""Apply Smith-Waterman on two AA sequences (local alignment).
	Using linear gap model.
	"""
	
	#Verify the arguments
	if not areArgumentsValid(seqA,seqB,gap,match,mismatch,matrixName):
		return
	gap=int(gap)
	if match!=None and mismatch!=None:
		match = int(match)
		mismatch = int(mismatch)
	
	print "Aligning sequence "+seqA+" with sequence "+seqB+" using Smith-Waterman:"
		
	#Compute the matrix of pair scores
	if matrixName==None:
		pairScores = createMatrixOfPairScores(seqA,seqB,match,mismatch)
	else:
		blosum50 = readScoringMatrix(matrixName)
		pairScores = createMatrixOfPairScores(seqA,seqB,matrix=blosum50)
	print "Here is the matrix of pair scores for the two sequences: "
	pprint(pairScores)
	
	#Generate the dynamic programming matrix (local align., linear gap model)
	# - initialize the first row of the dp matrix with zeros
	dpMatrix = []
	dpMatrix.append([0]*(len(seqA)+1))
	dpMatrixTrack = []
	dpMatrixTrack.append([trace[2]]*(len(seqA)+1)) # put '-' on first row
	overallMax = [-1,[]] # compute this maximum while generating the dp matrix 
	for i,aB in enumerate(seqB):
		dpMatrix.append([])
		dpMatrix[i+1].append(0) # - put zeros on first column
		dpMatrixTrack.append([])
		dpMatrixTrack[i+1].append(trace[1]) # - put '|' on first column
		for j,aA in enumerate(seqA):
			maxArgs = [dpMatrix[i][j] + pairScores[i][j],\
				dpMatrix[i][j+1] - gap,\
				dpMatrix[i+1][j] - gap,\
				0]
			maxVal = max(maxArgs)
			dpMatrix[i+1].append(maxVal)
			dpMatrixTrack[i+1].append(trace[maxArgs.index(maxVal)])
			if maxVal > overallMax[0]:
				overallMax[0] = maxVal
				overallMax[1] = (i+1,j+1)
	print "Here is the local dynamic programming matrix for the two sequences: "
	pprint(dpMatrix)
	print "Here is the local dynamic programming traceback for the two sequences: "
	pprint(dpMatrixTrack)
	
	#Backtracking
	# Start from the maximum value from the dp matrix
	loc = overallMax[1]
	alignment =[[],[]]
	# Go towards the top left corner until a zero value is found
	while dpMatrix[loc[0]][loc[1]] != 0:
		direction = dpMatrixTrack[loc[0]][loc[1]]
		if direction=='-':
			alignment[0].insert(0,seqA[loc[1]-1])
			alignment[1].insert(0,'-')
			loc = (loc[0],loc[1]-1)
		elif direction=='|':
			alignment[0].insert(0,'-')
			alignment[1].insert(0,seqB[loc[0]-1])
			loc = (loc[0]-1,loc[1])
		elif direction=='\\':
			alignment[0].insert(0,seqA[loc[1]-1])
			alignment[1].insert(0,seqB[loc[0]-1])
			loc = (loc[0]-1,loc[1]-1)		
	print "Here is a local alignment:"
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

if len(sys.argv)!=1 and len(sys.argv)!=6 and len(sys.argv)!=7:
	print "Usage: python align.py"
	print "Usage: python align.py global blossum50 gap HEAGAWGHEE PAWHEAE"
	print "Usage: python align.py global match mismatch gap HEAGAWGHEE PAWHEAE"
	exit(0)
elif len(sys.argv)==1:
	#Default behaviour
	seqA = "HEAGAWGHEE"
	seqB = "PAWHEAE"
	#needle(seqA,seqB,gap=2,match=1,mismatch=-2)
	needle(seqA,seqB,gap=8,matrixName="blosum50")
	seqA = "GHGKKVADALTN"
	seqB = "GHKRLLT"
	needle(seqA,seqB,gap=8,matrixName="blosum50")
	water(seqA,seqB,gap=8,matrixName="blosum50")
elif len(sys.argv)==6 or len(sys.argv)==7:
	#Alignment according to the received arguments
	# global or local
	#   using matrix or match/mismatch values
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