#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.structure import *

#########################################
def main():

	db, cursor = abca4_connect()

	for [id, protein] in hard_landing_search(cursor, "select id, protein from variants"):
		pattern = re.match('(\D+)(\d+)\D', protein)
		if not pattern:
			print(protein)
			continue
		aa  = pattern.group(1)
		pos = pattern.group(2)
		print(id, protein, aa, pos, find_region(int(pos)))
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
