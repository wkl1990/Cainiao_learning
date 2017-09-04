#Description: Fasta file exploration
#Environment: python3
#Data: 2017-09-04
#Author: WKL

#modules
import sys
sys.path.remove("~/usr/local/lib/python2.7/site-packages") #remove not available lib
import pysam
import os
#directory
os.chdir("./Cainiao_2")
#init, open and calculate
hg19 = pysam.FastaFile("UCSC_hg19_chrAll.fasta")
dir(hg19)
list(zip(hg19.references, hg19.lengths))
#output
for seqName in hg19.references:
	seq = hg19.fetch(seqName)
	A = seq.count("A") + seq.count("a")
	T = seq.count("T") + seq.count("t")
	C = seq.count("C") + seq.count("c")
	G = seq.count("G") + seq.count("g")
	N = seq.count("N") + seq.count("n")
	print(seqName, A, T, C, G, N)
