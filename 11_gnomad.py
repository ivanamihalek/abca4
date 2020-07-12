#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *



#########################################
def main():

	db, cursor = abca4_connect()

	for [variant_id, gdna_start, nt_from, nt_to, protein] in hard_landing_search(cursor, "select id, gdna_start, nt_from, nt_to, protein from variants"):
		if not gdna_start: continue
		qry  = "select position, reference, variant, variant_count, total_count, homozygote_count "
		qry += f"from gnomad.freqs_chr_1 where position={gdna_start}"
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue
		print(gdna_start, nt_from, nt_to, "      ", protein)
		for line in ret:
			[position, reference, variant, variant_count, total_count, homozygote_count] = line
			freq = "%.1e"%(float(variant_count)/total_count)
			print("\t", position, reference, variant, freq, homozygote_count)
			if reference==nt_from and variant==nt_to:
				qry = f"update variants set gnomad_freq={freq}, gnomad_homozygotes={homozygote_count} where id={variant_id}"
				ret = error_intolerant_search(cursor, qry)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
