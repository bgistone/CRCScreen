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
	--output	<file>	output file
	--help	call help
	'''
	exit()

opts, args = getopt.getopt(sys.argv[1:], 'hf:r:o', ['read1=', 'read2=', 'output='])
(read1, read2, output) = ("", "", "")
(total_reads, regular_reads, noPattern_reads, errorPattern_reads) = (0, 0, 0, 0)
(max, cross_matched, cross_unmatched) = (0, 0, 0)

for op, va in opts:
	if op in ['-h', '--help']:
		help()
	elif op in ['--read1', '-f']:
		read1 = va
	elif op in ['--read2', '-r']:
		read2 = va
	elif op in ['--output', '-o']:
		output = va


if not read1 and not read2:
	help()

if not output:
	output = "sample"

genotype={}
freq={}

shift = {"ACGTGCTG":1, "ATTGCTAC":2, "TATGTGCC":-1, "AACTGTGA":-2, "GGCACAT":1, "TCACAGT":2, "CAGCACG":-1, "GTGCGCAA":-2}
myDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
class Core:
	def __init__(self, np, ep, rr, npr, epr):
		self.NP = np
		self.EP = ep
		self.RR = rr
		self.NPR = npr
		self.EPR = epr
	def getGenotype_forSingle(self, ss):
		key = ""
		try:
			if read1:
				m = re.search(r'TGAATAT(\w)A(\w)\w{18}(\w{2})\w(\w{2}})\w(\w)\w+GATAC(\w)\w{35}(\w)',  ss)
			else:
				m = re.search(r'ATTGTT(\w)\w+AATTAGC(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{13}ACAAG(\w)',  ss)
			key = "".join(m.groups())
			key = self.__shift_check_forSingle(key, ss)
			self.RR += 1
		except AttributeError, e:
			try:
				if read1:
					m = re.search(r'AGG\w{22}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', ss)
				else:
					m = re.search(r'AAGAA\w{51}(\w)\w{35}(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{18}(\w)', ss)
				key = "".join(m.groups())
				key = self.__shift_check_forSingle(key, ss)
				self.RR += 1
			except AttributeError, e:
				self.NP.write("NP:\t" + ss+"\n")
				self.NPR += 1
		return key

	def getGenotype(self, s1, s2):
		key = ""
		try:
			m1 = re.search(r'TGAATAT(\w)A(\w)\w{18}(\w{2})\w(\w{2}})\w(\w)\w+GATAC(\w)\w{35}(\w)',  s1)
			m2 = re.search(r'ATTGTT(\w)\w+AATTAGC(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{13}ACAAG(\w)',  s2)
			k1 = "".join(m1.groups())
			k2 = "".join(m2.groups())
			#key = k1 + "/" + k2
			key = self.__shift_check(k1, k2, s1, s2)
			self.RR += 2
		except AttributeError, e:
			try:
				m1 = re.search(r'AGG\w{22}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', s1)
				m2 = re.search(r'AAGAA\w{51}(\w)\w{35}(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{18}(\w)', s2)
				k1 = "".join(m1.groups())
				k2 = "".join(m2.groups())
				#key = k1 + "/" + k2
				key = self.__shift_check(k1, k2, s1, s2)
				self.RR += 2
			except AttributeError, e:
				self.NP.write("NP:\t" + s1+"\n")
				self.NP.write("NP:\t" + s2+"\n")
				self.NPR += 2
		return key

	def __shift_check_forSingle(self, k, s):
		m = ""
		if k in shift:
			if read1:
				if shift.get(k) == 1:
					m = re.search(r'AGG\w{21}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', s)
				elif shift.get(k) == 2:
					m = re.search(r'AGG\w{20}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', s)
				elif shift.get(k) == -1:
					m = re.search(r'AGG\w{23}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', s)
				elif shift.get(k) == -2:
					m = re.search(r'AGG\w{24}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', s)
				try:
					k = "".join(m.groups())
				except AttributeError, e:
					self.EP.write("EP-read1:\t" + s + "\n")
					self.EPR += 1
		if read2:
			if k in shift:
				if shift.get(k) == 1:
					m = re.search(r'AAGAA\w{50}(\w)\w{35}(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{18}(\w)', s)
				elif shift.get(k) == 2:
					m = re.search(r'AAGAA\w{49}(\w)\w{35}(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{18}(\w)', s)
				try:
					k = "".join(m.groups())
				except AttributeError, e:
					self.EP.write("EP-read2:\t" + s + "\n")
					self.EPR += 1
		return k

	def __shift_check(self, k1, k2, s1, s2):
		m1 = ""
		m2 = ""
		if k1 in shift:
			#print k1 + "exist in shift"+shift.get(k1)+"\n"
			if shift.get(k1) == 1:
				m1 = re.search(r'AGG\w{21}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)', s1)
			elif shift.get(k1) == 2:
				m1 = re.search(r'AGG\w{20}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)',  s1)
			elif shift.get(k1) == -1:
				m1 = re.search(r'AGG\w{23}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)',  s1)
			elif shift.get(k1) == -2:
				m1 = re.search(r'AGG\w{24}(\w)A(\w)\w{18}(\w{2})\w(\w{2})\w(\w)\w{24}(\w)\w{35}(\w)',  s1)
			try:
				k1 = "".join(m1.groups())
			except AttributeError, e:
				self.EP.write("EP-read1:\t" + s1 + "\n")
				self.EPR += 1

		if k2 in shift:
			if shift.get(k2) == 1:
				m2 = re.search(r'AAGAA\w{50}(\w)\w{35}(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{18}(\w)', s2)
			elif shift.get(k2) == 2:
				m2 = re.search(r'AAGAA\w{49}(\w)\w{35}(\w)\w{24}(\w)\w(\w{2})\w(\w{2})\w{18}(\w)', s2)
			try:
				k2 = "".join(m2.groups())
			except AttributeError, e:
				self.EP.write("EP-read2:\t" + s2 + "\n")
				self.EPR += 1
		return k1 + "/" + k2

	def cross_match(self, str):
		myList = str.split('/')
		#i = len(myList[0]) - 1
		i = 1
		j = len(myList[1]) - 1
		while(j>=0):
			if myList[1][j] != myDict[myList[0][i]]:
				return 0			
			j = j - 1
			i = i + 1
		return 1

if __name__ == '__main__':
	(np, ep) = ("", "")
	if read1 and read2:
		np = open(output+".pair.unRegularReads", 'w')
		ep = open(output+".pair.errorPatternReads", 'w')
	elif read1:
		np = open(output+".read1.unRegularReads", 'w')
		ep = open(output+".read2.errorPatternReads", 'w')
	elif read2:
		np = open(output+".read2.unRegularReads", 'w')
		ep = open(output+".read2.errorPatternReads", 'w')

	toolBox = Core(np, ep, regular_reads, noPattern_reads, errorPattern_reads)
	if read1 and read2:
		print "ready go...."
		with open(read1, 'r') as finA, \
				open(read2, 'r') as finB:
					print "deal with: "+read1 +"\ndeal with: "+read2+"\n"
					for r1, r2 in izip(SeqIO.parse(finA, "fastq"), SeqIO.parse(finB, "fastq")):
						total_reads += 2
						seqR1 = str(r1.seq)
						seqR2 = str(r2.seq)
						key = toolBox.getGenotype(seqR1, seqR2)
						if key != "":
							genotype[key] = key
							if key not in freq:
								freq[key] = 1
							else:
								freq[key] += 1
		print "loading finished and ready to write result...\n"
		with open(output+".pair.CM.GT.freq", 'w') as fout1, \
				open(output+".pair.UM.GT.freq", 'w') as fout2:
					for k in freq:
						if max < freq[k]:
							max = freq[k]
						#str = str(genotype[k])
						r = toolBox.cross_match(genotype[k])
						string = str(genotype[k])+"\t"+str(freq[k])+"\n"
						if r == 1:
							if freq[k] >= 100:
								cross_matched += freq[k]
							fout1.write(string)
						else:
							fout2.write(string)
							if freq[k] >= 100:
								cross_unmatched+=freq[k]
	else:
		print "ready go...."
		if read1:
			read = read1
		else:
			read = read2
		with open(read, 'r') as fin:
			print "deal with: "+read+"\n"
			for r in (SeqIO.parse(fin, "fastq")):
				total_reads += 1
				seq = str(r.seq)
				key = toolBox.getGenotype_forSingle(seq)
				if key != "":
					genotype[key] = key
					if key not in freq:
						freq[key] = 1
					else:
						freq[key] += 1
		print "loading finished and ready to write result...\n"
		fout = ""
		if read1:
			fout = open(output+".read1.GT.freq",  'w')
		else:
			fout = open(output+".read2.GT.freq",  'w')
		for k in freq:
			if max < freq[k]:
				max = freq[k]
			if freq[k]>=100:
				cross_matched += freq[k]
			string = str(genotype[k]) + "\t" + str(freq[k]) + "\n"
			fout.write(string)
		fout.close()
	np.close()
	ep.close()

	#=================output log============
	print "\nReads loaded in total: ", total_reads
	print "no pattern found on: ", toolBox.NPR,", ", toolBox.EPR, " of these reads were detected frameshift events and can't been corrected\n"
	if read1 and read2:
		all = cross_matched + cross_unmatched
		mutation_count = cross_matched - max
		print cross_matched, "(", float(cross_matched)/all*100, "%) of ", all , " record were cross validated"
		print mutation_count, "(", float(mutation_count)/all*100, "%) mutation of ", all, "records were detected"
	else:
		all = cross_matched
		mutation_count = cross_matched - max
		print mutation_count, "(", float(mutation_count)/all*100,  "%) mutation of ", all, "records were detected"
		
	print "=====================All jobs done!===================\n"
