#!/usr/bin/env python
# coding=utf-8
import re, sys
import getopt
from itertools import *
from Bio import SeqIO

def help():
	print '''
	python genotype.py
	--read1	<file>	input read1
	--read2	<file>	input read2
	--help	call help
	'''
	exit()

opts, args = getopt.getopt(sys.argv[1:], 'hf:r', ['read1=', 'read2='])
(read1, read2) = ("", "")

for op, va in opts:
	if op in ['-h', '--help']:
		help()
	elif op in ['--read1', '-f']:
		read1 = va
	elif op in ['--read2', '-r']:
		read2 = va


if not read1 and not read2:
	help()



genotype={}
freq={}



with open(read1, 'r') as finA,\
		open(read2, 'r') as finB:
			for r1, r2 in izip(SeqIO.parse(finA, "fastq"), SeqIO.parse(finB, "fastq")):
				seqR1 = str(r1.seq)
				seqR2 = str(r2.seq)
				m1 = re.search(r'AGG\w{22}(\w)\w(\w)\w{14}(\w)\w{3}(\w{7})\w{24}(\w)', seqR1)
				m2 = re.search(r'CAC\w{9}\w{49}(\w)\w{35}(\w)\w{24}(\w{7})\w{3}(\w)\w{14}(\w)', seqR2)
				try:
					key = "".join(m1.groups())
					key = key + "/" +"".join(m2.groups())
					genotype[key] = key
					if key not in freq:
						freq[key] = 1
					else:
						freq[key] += 1
				except AttributeError, e:
					print seqR1, e
					print seqR2, e

with open(read1+".stat", 'w') as fout:
	for k in freq:
		string = str(genotype[k])+"\t"+str(freq[k])+"\n"
		fout.write(string)
