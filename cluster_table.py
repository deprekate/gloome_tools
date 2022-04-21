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

def css():
	return """
table { 
	/*	table-layout: auto; */
		table-layout: fixed;
        margin-left: 0em;
        margin-right: 0em;
        padding:1em 1em 1em 1em;
	    margin:1em 1em 1em 1em; 
        border-collapse: collapse;
      }
td {
        font-family: "Courier New", Courier, monospace;
        font-size:1em;
        font-weight: bold;
        text-align: center;
        overflow:hidden;
     /*   white-space:nowrap; */
		 white-space:pre;
/*	 	padding:0.1em 0.1em 0.1em 0.1em; */
/*	 	margin:0.5em 0.5em 0.5em 0.5em;*/
		width: 1em;
      }
td.Seq_Name{
        text-align: left;
        width: 15em;
        padding-right:1em;
      }
td.Score7{
		color: #FFFFFF;
        background: #B10026;
        }
td.Score6{
        background: #E31A1C;
        }
td.Score5{
        background: #FC4E2A;
        }
td.Score4{
        background: #FD8D3C;
	}
td.Score3{
        background: #FEB24C;
	}
td.Score2{
        background: #FED976;
        }
td.Score1{
        background: #B5B5B5;
        }
td.white{
	background: #FFFFFF;
"""

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
		idx2 = [i+1 for i in Z2['leaves'] ]


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
							if cell:
								f.write("\n")
					f.write("</tr>\n")
					A = None
					flag = False
					# add stuff below the bars
					'''
					f.write("<tr>\n")
					for cell in row:
						match = re.findall("title=\"[^\"]*", cell)
						if match:
							f.write("<td class='rotate'>\n")
							f.write(match[0].split("\"")[-1].split(":")[0])
						else:
							f.write("<td class='Seq_Name' style = 'text-align: right'>\n")
							f.write("og#")
						f.write("\n")
						f.write("</td>\n")
					f.write("</tr>\n")
					'''
				elif line.startswith('</table>'):
					f.write(line)
					flag = True
				elif line.startswith('<link rel='):
					pass
				elif flag:
					f.write(line)
				# This adds the css to keep the first column locked in place
				if line.startswith('</head>'):
					f.write('<style>\n')
					f.write('.Seq_Name {\n')
					f.write('  position: sticky;\n')
					f.write('  left: 0;\n')
					f.write('  background-color: #ffffff;\n')
					f.write('}\n')
					f.write(css())
					f.write('</style>\n')
					f.write(line)



