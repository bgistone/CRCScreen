#!/usr/bin/env python
# coding=utf-8
import re,  sys
import getopt
from collections import defaultdict
from itertools import *
from Bio import SeqIO

def help():
	print '''
	python genotype.py <genotype.list> <output>
	'''
	exit()

if len(sys.argv) < 3:
	help()

dict=['total', 'normal', '12S', '12R', '12C', '12D', '12A', '12V', '13R', '13C', '13D', '13V', '14I', 'nocountable']

all=defaultdict(defaultdict(dict))
all={}
info={}
samples=[]

with open(sys.argv[1], 'r') as LL:
	for line in LL:
		v = line.split('\t')
		sampleID=v[0]
		file = v[1].strip()
		samples.append(sampleID)

		with open(file, 'r') as IN:
			for data in IN.readlines():
				data.strip()
				arr = data.split()
				all.update({arr[0]:{sampleID:{'COUNT':arr[1], 'RATIO':arr[2]}}})

				info[arr[0]]=arr[3]
				print len(all)
		print "------------------", sampleID, "----------------"

with open(sys.argv[2], 'w') as OO:
	#title for output
	OO.write("GT\t")
	for s in samples:
		OO.write(s+"_freq\t"+s+"_ratio\t")
	OO.write("info\n")
	#data 
	for k in dict:
		OO.write(k+"\t")
		for sid in samples:
			if all.get(k):
				if all[k].get(sid):
					count = all[k][sid]['COUNT']
					ratio = all[k][sid]['RATIO']
					OO.write(count+"\t"+ratio+"\t")
				else:
					OO.write("0\t0\t")
			else:
				OO.write("0\t0\t")
		if info.get(k):
			OO.write(info[k]+"\n")
		else:
			OO.write("-\n")

