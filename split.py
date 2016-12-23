#!/usr/bin/env python
import sys
import getopt
import re
import gzip
from itertools import izip
from Bio import SeqIO
#from fuzzywuzzy import fuzz

PE = 0
SE = 0
total_read = 0

barcode = 'ATCTCG'
# function
def compile(barcode):
	pattern = '.*'+barcode+'.*'
	regex = re.compile(pattern)
	return regex
def help():
	print '''
	python split.py 
	--file1	<file>	input read1
	--file2	<file>	input read2
	--output	<file>	output file prefix
	--help		call for help
	'''
	exit()

regex = compile(barcode)
# parameters
opts, args = getopt.getopt(sys.argv[1:], 'hf:r:o', ['file1=', 'file2=', 'output='])
(file1, file2, output) = ("", "", "")

for op, va in opts:
	if op in ['-h', '--help']:
		help()
	elif op in ['--file1', '-f']:
		file1 = va
	elif op in ['--file2', '-r']:
		file2 = va
	elif op in ['--output', '-o']:
		output = va

if not file1 and not file2:
	help()
if not output:
	help()


with gzip.open(file1, 'r') as ff, \
		gzip.open(file2, 'r') as fr,\
		open(output+'_1.fq', 'w') as out_read1, \
		open(output+'_2.fq', 'w') as out_read2, \
		open('single.fq', 'w') as out_single:
	for fline, rline in izip(SeqIO.parse(ff, "fastq"), SeqIO.parse(fr, "fastq")):
		total_read += 2
		id1 = fline.description.split(':')
		id2 = rline.description.split(':')
		n1 = regex.search(id1[-1])
		n2 = regex.search(id2[-1])
		if n1 and n2:
			#out_read1.writelines([r1_id, r1_rd, r1_sn, r1_qa])
			#out_read2.writelines([r2_id, r2_rd, r2_sn, r2_qa])
			SeqIO.write(fline, out_read1, 'fastq')
			SeqIO.write(rline, out_read2, 'fastq')
			PE += 2
		elif n1:
			SE += 1
			#out_single.writelines([r1_id,  r1_rd,  r1_sn,  r1_qa])
			SeqIO.write(fline, out_single, 'fastq')
		elif n2:
			SE += 1
			#out_single.writelines([r1_id,   r1_rd,   r1_sn,   r1_qa])
			SeqIO.write(rline, out_single, 'fastq')

		if total_read % 10000 == 0:
			print "screened %d reads\n" % total_read

print "done with all reads\n"

with open("barcode_match_stat.xls", 'w') as OO:
	OO.writelines('barcode:%s\n' % barcode)
	OO.writelines('Pair End Reads:%d\n' % PE)
	OO.writelines('Single End Reads:%d\n' % SE)
	OO.writelines('Total Reads:%d\n' % total_read)
