#Description: GTF file exploration
#Environment: Python3
#Data: 2017-09-05
#Author: WKL

#module
import sys
import re
import os
import gzip
from collections import OrderedDict
from operator import itemgetter

#director 
os.chdir("./Cainiao_3")

#options, args
args = sys.argv

#init, definition
class Genome:
	"""docstring for Genome"""
	def __init__(self):
		self.chr = ""
		self.start = 0
		self.end = 0
class Gene(Genome):
	"""docstring for Gene"""
	def __init__(self):
		Genome.__init__(self)
		self.orientation = ""
		self.id = ""
class Transcript(Genome):
	"""docstring for Transcript"""
	def __init__(self):
		Genome.__init__(self)
		self.id = ""
		self.parent = ""
class Exon(Genome):
	"""docstring for Exon"""
	def __init__(self):
		Genome.__init__(self)
		self.parent = ""

def main(args):
	'''
	one input option
	input: a gtf file

	:return:

	'''
	chr_list = []
	gene_list = OrderedDict()
	transcript_list = OrderedDict()
	exon_list = []
	with gzip.open(args[1], 'rt') as fp_gtf:
		for line in fp_gtf:
			if line.startswith("#"):
				continue
			line = line.rstrip().split("\t")
			chr = line[0]
			type = line[2]
			start = int(line[3])
			end = int(line[4])
			orientation = line[6]
			attr = line[8]
			if not re.search(r'protein_coding', attr):
				continue
			if not chr in chr_list:
				chr_list.append(chr)
			if type == "gene":
				gene = Gene()
				id = re.search(r'gene_id "([^;]+)";?', attr).group(1)
				gene.chr = chr
				gene.start = start
				gene.end = end
				gene.id = id
				gene.orientation = orientation
				gene_list[id] = gene 
			elif type == "transcript":
				transcript = Transcript()
				id = re.search(r'transcript_id "([^;]+)";?', attr).group(1)
				parent = re.search(r'gene_id "([^;]+)";?', attr).group(1)
				if not parent in gene_list:
					continue
				transcript.chr = chr
				transcript.start = start 
				transcript.end = end 
				transcript.id = id 
				transcript.parent = parent 
				transcript_list[id] = transcript 
			elif type == "exon":
				exon = Exon()
				parent = re.search(r'transcript_id "([^;]+)";?', attr).group(1)
				if not parent in transcript_list:
					continue
				exon.chr = chr 
				exon.start = start 
				exon.end = end 
				exon.parent = parent 
				exon_list.append(exon)
	'''
	chr_gene(gene_list)
	gene_len(gene_list)
	gene_transcript(transcript_list)
	transcript_exon(exon_list)
	exon_pos(exon_list)
	'''
	gene_exon_pos(gene_list, transcript_list, exon_list)

def chr_gene(gene_list):
	'''
	gene number distribution in chromsome

	:param gene_list:
	:return:
	'''
	print("gene distribution in chromsome")
	count_gene = OrderedDict()
	for info in gene_list.values():
		chr = info.chr 
		if chr in count_gene:
			count_gene[chr] += 1
		else:
			count_gene[chr] = 1
	with open("chr_gene1.txt", 'w') as fp_out:
		for chr, num in count_gene.items():
			fp_out.write("\t".join([chr, str(num)]) + "\n")
			print("\t".join([chr, str(num)]) + "\n")

def gene_len(gene_list):
	'''
	gene length distribution

	:param gene_list:
	:return:
	'''
	print("gene length distribution")
	with open("gene_len1.txt", 'w') as fp_out:
		for gene_id, info in gene_list.items():
			len = info.end - info.start + 1
			fp_out.write("\t".join([gene_id, str(len)]) + "\n")
			print("\t".join([gene_id, str(len)]) + "\n")

def gene_transcript(transcript_list):
	'''
	transcript number distribution

	:param transcript_list:
	:return:
	'''
	print("transcript number distribution")
	count_tanscript = OrderedDict()
	for info in transcript_list.values():
		gene_id = info.parent
		if gene_id in count_tanscript:
			count_tanscript[gene_id] += 1
		else:
			count_tanscript[gene_id] = 1
	with open("gene_transcript1.txt", 'w') as fp_out:
		for gene_id, num in count_tanscript.items():
			fp_out.write("\t".join([gene_id, str(num)]) + "\n")
			print("\t".join([gene_id, str(num)]) + "\n")

def transcript_exon(exon_list):
	'''
	exon number statistics

	:param exon_list:
	:return:
	'''
	print("exon number distribution of transcript")
	count_exon = OrderedDict()
	for info in exon_list:
		transcript_id = info.parent
		if transcript_id in count_exon:
			count_exon[transcript_id] += 1
		else:
			count_exon[transcript_id] = 1
	with open("transcript_exon1.txt", 'w') as fp_out:
		for transcript_id, num in count_exon.items():
			fp_out.write("\t".join([transcript_id, str(num)]) + "\n")
			print("\t".join([transcript_id, str(num)]) + "\n")

def exon_pos(exon_list):
	'''
	exon position statistics
	
	:param exon_list:
	:return:
	'''
	print("exon position")
	count_exon = OrderedDict()
	for exon in exon_list:
		transcript_id = exon.parent
		start = exon.start
		end = exon.end 
		if transcript_id in count_exon:
			count_exon[transcript_id] += ",%s-%s" % (str(start), str(end))
		else:
			count_exon[transcript_id] = "%s-%s" % (str(start), str(end))
	with open("exon_pos.txt1", 'w') as fp_out:
		for transcript_id, pos in count_exon.items():
			fp_out.write("\t".join([transcript_id, pos]) + "\n")
			print("\t".join([transcript_id, pos]) + "\n")

def gene_exon_pos(gene_list, transcript_list, exon_list):
	'''
	exon 2 transcript 
	transcript 2 gene 
	gene 2 chromsome

	:param gene_list:
	:param transcript_list:
	:param exon_list:
	:return:
	'''
	print("gene 2 transcript 2 exon position")
	count_exon = OrderedDict()
	for exon in exon_list:
		transcript_id = exon.parent 
		start = exon.start 
		end = exon.end
		if transcript_id in count_exon:
			count_exon[transcript_id] += ", %s-%s" % (str(start), str(end))
		else:
			count_exon[transcript_id] = "%s-%s" % (str(start), str(end))
	count_transcript = OrderedDict()
	for transcript_id, info in transcript_list.items():
		gene_id = info.parent
		chr = info.chr 
		if gene_id not in count_transcript:
			count_transcript[gene_id] = OrderedDict()
		count_transcript[gene_id][transcript_id] = count_exon[transcript_id] 
	count_gene = OrderedDict()
	for gene_id, info in gene_list.items():
		chr = info.chr 
		if chr not in count_gene:
			count_gene[chr] = OrderedDict()
		count_gene[chr][gene_id] = count_transcript[gene_id]
	with open("gene_exon_pos.txt", 'w') as fp_out:
		for chr, gene in count_gene.items():
			for gene_id, transcript in gene.items():
				for transcript_id, exon in transcript.items():
					fp_out.write("\t".join([chr, gene_id, transcript_id, exon]) + "\n")

if __name__ == '__main__':	
	main(args)
		
		
