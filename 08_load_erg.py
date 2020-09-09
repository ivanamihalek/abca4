#!/usr/bin/python3

import os

from utils.mysql import *
from sys import argv


def parse_in(in_tsv):
	erg_string = {}
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()
	if not in_tsv[-4:] == ".tsv":
		print(f"is {in_tsv} a tsv file?")
		exit()
	inf = open(in_tsv)
	erg_vals = {}
	for line in inf:
		fields = line.strip().split('\t')
		if len(fields)<4: continue
		if 'Pending' in fields[2]: continue
		if not fields[0] or len(fields[0])==0: continue
		try:
			patient_id = int(fields[0])
		except:
			continue
		if not fields[3] or len(fields[3])==0: continue
		if patient_id not in erg_vals: erg_vals[patient_id] = []
		erg_vals[patient_id].append(f"{fields[1]}:{fields[3]}")

	for pid, vals in erg_vals.items():
		erg_string[pid] = ";".join(erg_vals[pid])

	return erg_string

#########################################
def main():

	db, cursor = abca4_connect()

	if len(argv) < 2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()
	tsv = argv[1]
	ergs = parse_in( tsv)
	for pid, ergstr in ergs.items():
		print(pid, ergstr)
		store_or_update(cursor, 'cases',
		                fixed_fields={'publication_id':52, "patient_xref_id":pid},
		                update_fields={'erg_r_rod':ergstr})
	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
