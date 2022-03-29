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

sys.path.pop(0)
from genbank.file import File

def is_valid_file(x):
	if x and not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def just_seq(text):
	seq = ''
	for line in text.split("\n"):
		if line and not line.startswith(">"):
			seq += line
	return seq
	

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('MICROBIALIZER_ZIP', type=is_valid_file, help='input file')
	parser.add_argument('GLOOME_DIR', type=is_valid_file, help='input file')
	args = parser.parse_args()


	# THIS IS TO GET THE INFORMATION FROM THE MICROBIALIZER ZIP
	orthologs = dict()
	with ZipFile(args.MICROBIALIZER_ZIP) as zf:
		with zf.open('11_final_table/final_orthologs_table.csv', 'r') as csv_file:
			csv_reader = csv.reader(TextIOWrapper(csv_file, 'utf-8'))
			csv_columns = next(csv_reader)
			for csv_row in csv_reader:
				# process the CSV here
				with zf.open('12_orthologs_groups_dna_sequences/'+csv_row[0]+'_dna.fas', 'r') as fas_file:
					orthologs[csv_row[0]] = ''
					for fas_row in fas_file:
						if fas_row.decode().startswith('>') and orthologs[csv_row[0]]:
							break
						orthologs[csv_row[0]] += fas_row.decode()

	# THI IS TO ANNOTATE THE GLOOME GAIN/LOSS HTML PAGES
	text = []
	i = 0
	with open(os.path.join(args.GLOOME_DIR, 'MSA_color_coded_by_gain_probability.html'), "r+") as f:
		for line in f:
			if line.startswith('<td valign = bottom'):
				line = line.replace('title=', 'oncontextmenu="navigator.clipboard.writeText(\'' + orthologs['og_'+str(i)].replace("\n", "\\n") + '\');return false;" title=')
				line = line.replace('title="', 'title="' + just_seq(orthologs['og_'+str(i)]) + ' right click to copy sequence ')
				i += 1
			text.append(line)

	with open(os.path.join(args.GLOOME_DIR, 'MSA_color_coded_by_gain_probability.html'), "w+") as f:
		for line in text:
			f.write(line)
