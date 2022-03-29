#!/usr/bin/env python3
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from itertools import chain as chain
import csv
from io import TextIOWrapper
from zipfile import ZipFile
import re

sys.setrecursionlimit(10000)

import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform

sys.path.pop(0)
from genbank.file import File

def is_valid_file(x):
	if x and not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def strip_html(text):
	try:
		return int(re.sub('<[^<]+?>', '', text))
	except:
		return 0

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('GLOOME_DIR', type=is_valid_file, help='input file')
	args = parser.parse_args()



	# THI IS TO ANNOTATE THE GLOOME GAIN/LOSS HTML PAGES
	for kind in ['gain','loss']:
		text = []
		i = 0
		flag = False
		table = []
		width = None
		with open(os.path.join(args.GLOOME_DIR, 'MSA_color_coded_by_'+kind+'_probability.html'), "r+") as f:
			for line in f:
				if line.startswith('</table>'):
					flag = True
				elif flag:
					pass
				elif line.startswith('<tr>'):
					if table and not width:
						width = len(table[-1])
					elif table and len(table[-1]) != width:
						table[-1] = [''] * width
					table.append([])
				elif line.startswith('<td '):
					table[-1].append(line.rstrip())
				text.append(line)
		
		fig = plt.figure(figsize=(8,8))
	
		if len(table[-1]) != width:
			table[-1] = [''] * width
		table = np.array(table)
	
		table_vec = np.vectorize(strip_html)
		table_num = table_vec(table)
		taxa = table_num.shape[0] - 5
		d2 = sch.distance.pdist(table_num[:taxa,1:].transpose())
		D2 = squareform(d2)
		Y2 = sch.linkage(d2, method='single')
		Z2 = sch.dendrogram(Y2)
		idx2 = Z2['leaves']
	
		d1 = sch.distance.pdist(table_num[:taxa,])
		D1 = squareform(d1)
		Y1 = sch.linkage(d1, method='single')
		Z1 = sch.dendrogram(Y1)
		idx1 = Z1['leaves']
	
		A = table[idx1+list(range(taxa,taxa+5)),:]
		A = A[:,[0]+idx2]	
	
		flag = True
		with open(os.path.join(args.GLOOME_DIR, 'MSA_color_coded_by_'+kind+'_probability.html'), "w+") as f:
			for line in text:
				if line.startswith('<table>') and A is not None:
					f.write(line)
					for row in A:
						f.write("<tr>\n")
						for cell in row:
							f.write(cell)
							f.write("\n")
					f.write("</tr>\n")
					A = None
					flag = False
				elif line.startswith('</table>'):
					flag = True
				elif flag:
					f.write(line)












