#!/usr/bin/env python
# coding=utf-8
import re,sys
import getopt
import gzip
#from fuzzysearch import find_near_matches
#from fuzzysearch import find_near_matches_with_ngrams
from fuzzywuzzy import fuzz
from itertools import *
from Bio import SeqIO

def help():
	print '''
	python genotype.py
	-f, --read1 <file>  input read1 fastq
	-r, --read2 <file>  input read2 fastq
	-o, --output    <file>  output file
	-c, --cutoff	<int>	cutoff of match length[20]
	--help  call help
	seed for fishing:
		read1: TAGGCAAGAGTGCCTTGACGATAC
		read2: GTATCGTCAAGGCACTCTTGCCTA
	'''
	exit()

opts,  args = getopt.getopt(sys.argv[1:],  'hf:r:o::c::',  ['read1=',  'read2=',  'output=', 'cutoff='])
(read1,  read2,  output, cutoff) = ("",  "",  "sample", 20)

for op,  va in opts:
	if op in ['-h',  '--help']:
		help()
	elif op in ['--read1',  '-f']:
		read1 = va
	elif op in ['--read2',  '-r']:
		read2 = va
	elif op in ['--output',  '-o']:
		output = va
	elif op in ['--cutoff', '-c']:
		cutoff = int(va)

if not read1 and not read2:
	help()

def NWDistance(seedSequence, candidateSequence):
	s=-1 #a mismatch would deduce 1 point.
	m=1 #plus 1 point for one match.
	g=-2 #deduce 2 point for one gap.  
	seedSequence=seedSequence.strip()
	candidateSequence=candidateSequence.strip()
	if len(seedSequence)==0:
		print "Error,  seed sequence length equal zero."
		sys.exit(1)
	elif len(candidateSequence)==0:
		print "Error,  candidate sequence length equal zero."
		sys.exit(1)
	sLen=len(seedSequence)
	cLen=len(candidateSequence)
	table=[]
	for m in range(0, len(seedSequence)+1):
		table.append([m*g])
	table[0]=[]
	for n in range(0, len(candidateSequence)+1):
		table[0].append(n*g)
	for i in range(sLen):
		for j in range(cLen):
			table[i+1].append(max(table[i][j]+(m if seedSequence[i] == candidateSequence[j] else s),table[i][j+1]+g, table[i+1][j]+g))	
	#	
	i=sLen-1
	j=cLen-1
	NewSeed = seedSequence[i]
	NewCandidate = candidateSequence[j]
	if len(seedSequence) <= 1 or len(candidateSequence)<=1:
		print "Error,  too short!"
		sys.exit(1)
	while True:					
		if i == 0 and j == 0:
			break
		if seedSequence[i] == candidateSequence[j]:
			if table[i-1][j-1]+1>table[i-1][j]-2 and table[i-1][j-1]+1>table[i][j-1]-2:
				i = i - 1
				j = j -1
				NewSeed = u"%s%s" % (seedSequence[i],  NewSeed)
				NewCandidate = u"%s%s" % (candidateSequence[j],  NewCandidate)
			else:
				if table[i][j+1] > table[i+1][j]:
					i = i-1
					NewSeed = u"%s%s" % (seedSequence[i],  NewSeed)
					NewCandidate = u"%s%s" % ('-',  NewCandidate)
				else:
					j = j-1
					NewSeed = u"%s%s" % ('-',  NewSeed)
					NewCandidate = u"%s%s" % (candidateSequence[j],  NewCandidate)
		else:
			if table[i-1][j-1]+1>table[i-1][j]-2 and table[i-1][j-1]+1>table[i][j-1]-2:
				i = i - 1
				j = j - 1
				NewSeed = u"%s%s" % (seedSequence[i],  NewSeed)
				NewCandidate = u"%s%s" % (candidateSequence[j],  NewCandidate)
			else:
				if table[i][j+1] > table[i+1][j]:
					i = i-1
					NewSeed = u"%s%s" % (seedSequence[i],  NewSeed)
					NewCandidate = u"%s%s" % ('-',  NewCandidate)
				else:
					j = j-1
					NewSeed = u"%s%s" % ('-',  NewSeed)
					NewCandidate = u"%s%s" % (candidateSequence[j],  NewCandidate)
	#distance
	mismath=0
	math=0
	gap=0
	charZipList=zip(NewSeed, NewCandidate)
	#delete the head gap
	for n in range(len(charZipList)):
		if "-" in charZipList[n]:
			del charZipList[0]
		else:
			break
	#delete the tail gap
	while True:
		lastTuple=charZipList.pop()
		if "-" in lastTuple:
			continue
		else:
			charZipList.append(lastTuple)
			break
	#
	for n in range(len(charZipList)):
		charTuple=charZipList[n]
		if charTuple[0]==charTuple[1]:
			math+=1
		elif "-" in charTuple:
			gapLoc=charTuple.index("-")
			if charZipList[n+1][gapLoc]=="-":
				continue
			else:
				gap+=1
		else:
			mismath+=1
	distance=round(1.0-float(math)/float(mismath+math+gap), 4)
	return distance

def align(seed, read, start, end):
	seq = read[start:end]
	slen = len(seed)
	clen = len(seq)
	dis = clen - slen;
	max_match_len = 0;
	for i in range(dis):
		match_len = 0
		for j in range(slen):
			if seed[j] == seq[i+j]:
				match_len += 1
		if match_len > max_match_len:
			max_match_len = match_len
	return max_match_len

if __name__ == '__main__':
	subseq_for_read1 = "TAGGCAAGAGTGCCTTGACGATAC"
	subseq_for_read2 = "GTATCGTCAAGGCACTCTTGCCTA"
	methy_subseq_for_read1 = "TAGGCAAGAGTGCCTTGACGATAC"
	methy_subseq_for_read2 = "GTATCGTCAAGGCACTCTTGCCTA"
	with gzip.open(read1, 'rb') as finA,\
			gzip.open(read2, 'rb') as finB,\
				open(output+"_filter_1.fastq", 'w') as foutA,\
					open(output+"_filter_2.fastq", 'w') as foutB:
						for r1, r2 in izip(SeqIO.parse(finA, 'fastq'), SeqIO.parse(finB, 'fastq')):
							#fm1=find_near_matches(subseq_for_read1, str(r1.seq), max_l_dist=1)
							#fm2=find_near_matches(subseq_for_read2, str(r2.seq), max_l_dist=1)
							#fm1 = fuzz.ratio(str(r1.seq), subseq_for_read1)
							#fm2 = fuzz.ratio(str(r2.seq), subseq_for_read2)
							#fm1 = NWDistance(subseq_for_read1, str(r1.seq))
							#fm2 = NWDistance(subseq_for_read2, str(r2.seq))
							fm1 = align(subseq_for_read1, r1.seq, 55, 87)
							fm2 = align(subseq_for_read2, r2.seq, 96, 128)
							#print str(r1.seq)+"\tTAGGCAAGAGTGCCTTGACGATAC\t"+str(fm1)
							#print str(r2.seq)+"\tGTATCGTCAAGGCACTCTTGCCTA\t"+str(fm2)
							#print "--", fm1, fm2, "--"
							if (fm1 >= cutoff) and (fm2 >= cutoff):
								SeqIO.write(r1, foutA, 'fastq')
								SeqIO.write(r2, foutB, 'fastq')
							else:
								fm1 = align(methy_subseq_for_read1,  r1.seq,  55,  87)
								fm2 = align(methy_subseq_for_read2,  r2.seq,  96,  128)
								if (fm1 >= cutoff) and (fm2 >= cutoff):
									SeqIO.write(r1,  foutA,  'fastq')
									SeqIO.write(r2,  foutB,  'fastq')
