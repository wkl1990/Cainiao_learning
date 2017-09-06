#Description: merge individual count files to one matrix
#Environment: python3
#Encoding: UTF-8
#Data: 2017-09-06
#Author: WKL

#modules
import re
import os
import gzip
import glob
import time
from collections import OrderedDict
#time start
start = time.clock()
#directory
path = "./Cainiao_4"
os.chdir(path)
#init, definition
count_matrix = OrderedDict()
def Readfile(File):
	if re.search(r'\.gz$', File):
		fp_in = gzip.open(File, 'rt')
	else:
		fp_in = open(File, 'r')
	try:
		while True:
			line = fp_in.readline()
			if not line:
				break
			line = line.rstrip()
			array = line.split("\t")
			gene_name = array[0]
			gene_count = array[1]
			if gene_name in count_matrix:
				count_matrix[gene_name].append(gene_count)
			else:
				count_matrix[gene_name] = [gene_count]
	finally:
		fp_in.close()
	return count_matrix
#main, output
def main():
	file_list = glob.glob("./count_files/*.txt*")
	for file in file_list:
		Readfile(file)
	with open("count_matrix2.txt", 'w') as fp_out:
		for gene_name, gene_counts in count_matrix.items():
			fp_out.write("%s\t%s\n" % (gene_name, "\t".join(gene_counts)))

if __name__ == '__main__':
	main()
#time end
end = time.clock()
print("used %s s" % str(end - start) )
