#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.structure import *

#########################################
def main():

	db, cursor = abca4_connect()

	for [variant_id, protein] in hard_landing_search(cursor, "select id, protein from variants"):
		pattern = re.match('(\D+)(\d+)\D', protein)
		if not pattern:
			print(protein)
			continue
		aa  = pattern.group(1)
		pos = pattern.group(2)
		structural_domain = find_region(int(pos))
		print(variant_id, protein, aa, pos, structural_domain)
		qry = f"update variants set protein_domain='{structural_domain}' where id={variant_id}"
		error_intolerant_search(cursor, qry)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
