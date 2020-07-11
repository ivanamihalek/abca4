#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *



#########################################
def main():
	# parse_cdna('972_973delinsAT')
	# # #print(get_cdna()[3810])
	# # #parse_cdna('4510_4535del')
	# # #print(get_cdna()[4509:4535])
	# exit()

	db, cursor = abca4_connect()

	for [gdna_start, nt_from, nt_to, protein] in hard_landing_search(cursor, "select gdna_start, nt_from, nt_to, protein from variants"):
		if not gdna_start: continue
		qry = f"select position, reference, variant, variant_count, total_count from gnomad.freqs_chr_1 where position={gdna_start}"
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue
		print(gdna_start, nt_from, nt_to, "      ", protein)
		for line in ret:
			[position, reference, variant, variant_count, total_count] = line
			print("\t", position, reference, variant, "%.1e"%(float(variant_count)/total_count))

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
