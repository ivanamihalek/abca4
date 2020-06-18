#!/usr/bin/python3

import os
from utils.mysql import *
from sys import argv


def parse_in(in_tsv):
	inf = open(in_tsv)
	for line in inf:
		fields = [f.strip() for f in line.split('\t')]
		if 'pubmed' in fields[0].lower(): continue # this is header
		[pmid, ref, patient_id, cdna1, protein1, cdna2, protein2, hapl_tested, comment, value_type, onset] = fields[:11]
		progression = [f for f in fields[11:]]
		if len(progression)==0:
			print("progression not given", pmid, ref, patient_id)
			exit()

		if len(progression)%3!=0:
			print("length of progression not divisible by 3", pmid, ref, patient_id)
			exit()

		if hapl_tested.lower() != 'yes':
			continue # let's stay away from those for now

		print([pmid, patient_id, cdna1, protein1, cdna2, protein2, value_type, onset])

	inf.close()

#########################################
def main():

	if len(argv)<2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()

	parse_in(in_tsv)

	db,cursor = abca4_connect()


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

