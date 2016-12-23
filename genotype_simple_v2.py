#!/usr/bin/env python
# coding=utf-8
import re, sys
import getopt
from itertools import *
from Bio import SeqIO

def help():
	print '''
	python genotype.py
	-f, --read1	<file>	input read1
	-r, --read2	<file>	input read2
	-o, --output	<file>	output file[sample]
	-c, --cutoff	<int>	cutoff of sub-genotype frequence for stat[100]
	-h, --help	call help
	'''
	exit()

opts, args = getopt.getopt(sys.argv[1:], 'hf:r:o:c:', ['read1=', 'read2=', 'output=', 'cutoff'])
(read1, read2, output, cutoff) = ("", "", "sample", 100)
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
	elif op in ['--cutoff', '-c']:
		cutoff = int(va)

if not read1 and not read2:
	help()

#================init====================
genotype={}
freq={}
disease_ratio = {}
anno = {}
shift = {
		"ACGTGCTG":1, 
		"ATTGCTAC":2, 
		"TATGTGCC":-1,
		"AACTGTGA":-2,
		"GGCACAT":1, 
		"TCACAGT":2, 
		"CAGCACG":-1, 
		"GTGCGCAA":-2
		}
myDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
#AAGGGGGAC
reference = {0:'A', 1:'A', 2:'G', 3:'G', 4:'G', 5:'G', 6:'G', 7:'A', 8:'C'}

code = {'2T':'12C', 
		'2C':'12R',
		'2A':'12S',
		'3C':'12A', 
		'3A':'12D', 
		'3T':'12V',  
		'4T':'13C', 
		'4C':'13R',
		'5A':'13D',
		'5T':'13V', 
		'6A':'14I'
		}

location = {0:'12', 1:'14', 2:'34', 3:'35', 4:'37', 5:'38', 6:'40', 7:'65', 8:'81'}

aa = {'12C':'(p.Gly12Cys)', 
		'12R':'(p.Gly12Arg)',
		'12S':'(p.Gly12Ser)',
		'12A':'(p.Gly12Ala)',
		'12D':'(p.Gly12Asp)', 
		'12V':'(p.Gly12Val)', 
		'13D':'(p.Gly13Asp)', 
		'13V':'(p.Gly13Val)', 
		'13C':'(p.Gly13Cys)', 
		'13R':'(p.Gly13Arg)', 
		'14I':'(p.Val14Ile)'}

#================core===============
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
		var = ""
		i = 1
		j = len(myList[1]) - 1
		while(j>=0):
			if myList[1][j] != myDict[myList[0][i]]:
				return 0			
			#reference[index] = loci
			#exon[index] = exon_loci
			#if myList[0][i] != reference[i]:
			#	var = var + "\t" + exon[i] + ":" + reference[i] + "->" + myList[0][i]
			j = j - 1
			i = i + 1
		#if not var.strip():
		#	return var
		return 1

	def compliment(self, key):
		newKey = ""
		for i in range(len(key)):
			newKey = myDict[key[i]] + newKey
		return "A"+newKey

	def annotation(self, key):
		myList = key.split('/')
		for i in range(2, 7):
			if myList[0][i] != reference[i]:
				tmp = str(i) + myList[0][i]
				
				if tmp in code:
					var_key = code[tmp]
					annotation = "NM_033360(KRAS):c." + location[i] + reference[i] + "->" + myList[0][i] + aa[var_key]
					anno[var_key] = annotation
					if var_key not in disease_ratio:
						disease_ratio[var_key] = freq[key]
					else:
						disease_ratio[var_key] += freq[key]
				#else:
				#	annotation = "NM_033360(KRAS):c." + location[i] + reference[i] +"->" + myList[0][i] + "(-)"
					
				

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
	log = open(output+".log", 'w')

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
		print "variation searching done\n\n"
		print "dump genotype and freqence:\n\t" + output + ".pair.CM.GT.freq\n\t" + output + ".pair.UM.GT.freq\n"
		with open(output+".pair.CM.GT.freq", 'w') as fout1, \
				open(output+".pair.UM.GT.freq", 'w') as fout2:
					for k in freq:
						if max < freq[k]:
							max = freq[k]
						#r = toolBox.cross_match(genotype[k])
						r = toolBox.cross_match(k)
						string = str(genotype[k])+"\t"+str(freq[k])+"\n"
						if r == 1:
							if freq[k] >= cutoff:
								cross_matched += freq[k]
								#variant annotation
								toolBox.annotation(k)
							fout1.write(string)
						else:
							fout2.write(string)
							if freq[k] >= cutoff:
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
	
		print "variants searching done\n\n"
		
		fout = ""
		if read1:
			fout = open(output+".read1.GT.freq",  'w')
			print "dump genotype and freqence:\n\t" + output + ".read1.GT.freq\n"
		else:
			fout = open(output+".read2.GT.freq",  'w')
			print "dump genotype and freqence:\n\t" + output + ".read2.GT.freq\n"
		for k in freq:
			if max < freq[k]:
				max = freq[k]
			if freq[k]>=cutoff:
				cross_matched += freq[k]
				toolBox.annotation(k)
			string = str(genotype[k]) + "\t" + str(freq[k]) + "\n"
			fout.write(string)

			#variant annotation
			if read1:
				toolBox.annotation[k]
			else:
				key = toolBox.compliment(k)
				toolBox.annotation[key]
		fout.close()
	np.close()
	ep.close()
	#=================variation=============
	with open(output+".disease_variants.anno.xls", "w") as fout:
		fout.write("total\t" + str(cross_matched) +"\t100"+"\ttotal\n")
		fout.write("normal\t" + str(max) + "\t" + str(float(max)/cross_matched * 100)+"\tnormal\n")
		uncountable = cross_matched - max;
		for vk in disease_ratio:
			uncountable -= disease_ratio[vk]
			ratio = float(disease_ratio[vk]) / cross_matched * 100
			fout.write(vk + "\t" + str(disease_ratio[vk]) + "\t" + str(ratio) + "\t" + anno[vk] + "\n")
		fout.write("uncountable\t"+str(uncountable)+"\t"+str(float(uncountable)/cross_matched*100)+"\tuncountable\n")
	#=================output log============
	log.write("Reads loaded in total: "+str(total_reads)+"\n")
	log.write("no pattern found on: "+str(toolBox.NPR)+", "+str(toolBox.EPR)+ " of these reads were detected frameshift events and can't been corrected\n")
	if read1 and read2:
		all = cross_matched + cross_unmatched
		mutation_count = cross_matched - max
		log.write(str(cross_matched)+ "("+ str(float(cross_matched)/all*100)+ "%) of "+ str(all)+ " record were cross validated\n")
		log.write(str(mutation_count)+ "("+ str(float(mutation_count)/all*100)+ "%) mutation of "+ str(all)+ " records were detected\n")
	else:
		all = cross_matched
		mutation_count = cross_matched - max
		log.write(str(mutation_count)+ "("+ str(float(mutation_count)/all*100)+  "%) mutation of "+ str(all)+ " records were detected\n")
	log.close()
	print "=====================All jobs done!===================\n"
