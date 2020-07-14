#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.structure import *

# TODO collect orthod and paras, run specs
# read in
# store to table - the variant columns are tinyints called  conserved_in_ortho and conserved_in_para
def read_specs(infile):
	if not os.path.exists(infile):
		print(infile, "not found")
		exit()

	cons = {}
	inf = open(infile, "r")
	header =  None
	for line in inf:
		if not header:
			# the first field is the comment sign
			header = line.strip().split()[1:]
		else:
			named_field = dict(zip(header,line.strip().split()))
			if float(named_field["rvet"])<0.23:
				cons[int(named_field["pos_in_GNAO1"])] = "yes"

	inf.close()

	return cons

#########################################
def main():

	db, cursor = abca4_connect()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
