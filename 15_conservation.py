#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.annotation import single_letter_code


# store to table - the variant columns are tinyints called  conserved_in_ortho and conserved_in_para
def read_specs(infile, cutoff):
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
			if float(named_field["rvet"])<cutoff:
				cons[int(named_field["pos_in_human"])] = named_field["human"]

	inf.close()

	return cons


def fill_column(cursor, filename, column, cutoff):
	cons = read_specs(filename, cutoff)

	for [variant_id, protein] in hard_landing_search(cursor, "select id, protein from variants"):
		pattern = re.match('(\D+)(\d+)\D', protein)
		if not pattern: continue
		aa  =  single_letter_code[pattern.group(1).upper()]
		pos = int(pattern.group(2))
		conserved = cons.get(pos, False)
		if not conserved: continue
		if conserved!=aa:
			print(variant_id, protein, aa, pos, conserved,  "mismatch")
			exit()
		print(variant_id, protein, aa, pos, cons.get(pos, "not conserved"))
		qry = f"update variants set {column}=1 where id={variant_id}"
		error_intolerant_search(cursor, qry)


#########################################
def main():

	db, cursor = abca4_connect()
	fill_column(cursor, "conservation/abca4/vertebrates_only/specs_out.score", "conserved_in_ortho_verts", 0.40)
	fill_column(cursor, "conservation/abca4/paras/specs_out.score", "conserved_in_para_verts", 0.16)

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
