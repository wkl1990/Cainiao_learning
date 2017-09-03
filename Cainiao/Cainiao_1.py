##Description: Calculate the length of human exons
##Envrionment: Python3
##Data: 2017-09-03
##Author: WKL

#import modules
import re
import os
from collections import OrderedDict
from operator import itemgetter
#change directory
os.chdir("./Cainiao_1")
#init
exonLength = 0
overlapExon = OrderedDict()
#open
with open("CCDS.20160908.txt",'rt') as f:
	for line in f:
		if line.startswith("#"):
			continue
		line = line.rstrip()
		lst = line.split("\t")
		if lst[-2] == "-":
			continue
		lst[-2] = re.sub("\[|\]","",lst[-2])
		exons = lst[-2].split(", ")
		for exon in exons:
			start = int(exon.split("-")[0])
			end = int(exon.split('-')[1])
			coordinate = lst[0]+":"+exon
			if coordinate not in overlapExon.keys():
				overlapExon[coordinate] = 1
				exonLength += end - start

#output 
print(exonLength)
