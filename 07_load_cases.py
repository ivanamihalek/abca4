#!/usr/bin/python3

import os
import re

from utils.mysql import *
from sys import argv, stdout

def parse_in(in_tsv):
	cases  = []
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()
	if not in_tsv[-4:]==".tsv":
		print(f"is {in_tsv} a tsv file?")
		exit()
	inf = open(in_tsv)
	for line in inf:
		fields = [f.strip() for f in line.split('\t') if f!="blah"]
		print(fields)
		cases.append(fields)

	return cases

def  store_variants(cursor, c, p):
	print(f"storing {c} {p}")
	cdna_vars = [v.replace(" ","") for v in c.split(";")]
	protein_vars = [v.replace(" ","") for v in p.split(";")]
	if len(cdna_vars)==0:
		print(f"cdna vars empty  {c} {p}")
		exit()
	if len(protein_vars)==0:
		print(f"proteins vars empty  {c} {p}")
		exit()

	if len(cdna_vars)!=len(protein_vars):
		print(f"var lengths not equal  {c} {p}")
		exit()

	for i in range(len(cdna_vars)):
		cdna_var = cdna_vars[i]
		protein_var = protein_vars[i]


#########################################
def main():

	if len(argv)<2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	cases = parse_in(in_tsv)

	db,cursor = abca4_connect()

	for case  in cases:
		[pubmed_id,  patient_id, c1, p1, c2, p2, value, onset, progression_string] = case
		# for allele1, allele 2
		allele_ids = []
		for [cdna_variants, protein_variants] in [[c1,p1], [c2,p2]]:
			# store or retrieve variant ids
			variant_ids = store_variants(cursor, c1, p1)
			# store or retrieve allele id
			#allele_ids.append(store_allele(cursor, variant_ids))
		# store or retreieve publication id
		#publication_id = store_publication(pubmed_id)
		# store case: allele_id_1, allele_id_2, publication_id, patient_id, onset, value, better, progression
		# store_or_update(cursor, fixed_fields={'publication_id': publication_id, 'patient_xref_id': patient_id},
		# 						update_fields ={'allele_id_1': allele_ids[0], 'allele_id_2': allele_ids[1],
		# 										'onset_age': onset, 'acuity_type':'BCVA', 'eye':'better', 'progression':progression_string})

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

