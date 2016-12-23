#!/usr/bin/env python
import re, sys

genotype={}
freq={}
with open("TAGGCAAGAGTGCCTTGACGATAC.r1.fq", 'r') as fin:
	for line in fin.readlines():
		#m = re.search(r'TC(\w)\w\w\wGTT(\w)GATCATA', line)
		#m1 = re.search(r'GC(\w).*TTGCCTA(\w{7})AGC', line)
		#m2 = re.search(r'(\w)CCAAC\wACCACAAG(\w)', line)
		line = fin.readlines()
		m = re.search(r'ACGAGG\w{22}(\w)\w(\w)\w{14}(\w)\w{3}(\w{7})\w{24}(\w)', line)
		#print m.groups()
		key = "".join(m.groups())
		genotype[key] = m.groups()
		if key not in freq:
			freq[key] = 1
		else:
			freq[key] += 1
		line = fin.readlines()
		line = fin.readlines()

with open("TAGGCAAGAGTGCCTTGACGATAC.r1.fq.stat", 'w') as fout:
	for k in freq:
		string = str(genotype[k])+"\t"+str(freq[k])+"\n"
		fout.write(string)

genotype={}
freq={}
with open("GTATCGTCAAGGCACTCTTGCCTA.r2.fq", 'r') as fin:
	for line in fin.readlines():
		#m = re.search(r'TC(\w)\w\w\wGTT(\w)GATCATA', line)
		#m1 = re.search(r'GC(\w).*TTGCCTA(\w{7})AGC', line)
		#m2 = re.search(r'(\w)CCAAC\wACCACAAG(\w)', line)
		line = fin.readlines()
		m = re.search(r'CAC\w{9}\w{49}(\w)\w{35}(\w)\w{24}(\w{7})\w{3}(\w)\w{14}(\w)', line)
		#print m.groups()
		key = "".join(m.groups())
		genotype[key] = m.groups()
		if key not in freq:
			freq[key] = 1
		else:
			freq[key] += 1
		line = fin.readlines()
		line = fin.readlines()

with open("GTATCGTCAAGGCACTCTTGCCTA.r2.fq.stat", 'w') as fout:
	for k in freq:
		string = str(genotype[k])+"\t"+str(freq[k])+"\n"
		fout.write(string)

