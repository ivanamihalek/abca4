#!/usr/bin/python3

import os
import re

from utils.mysql import *
from sys import argv, stdout
from utils.abca4_gene import *


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
		cases.append(fields)
	return cases


###########################################################
def gdna_from_cdna(cdna_pos):
	if cdna_pos is None:
		print("cdna cannot be none here (do something)")
		exit()
	offset = 0
	if "+" in cdna_pos:
		[cdna_pos, offset] = [int(i) for i in cdna_pos.split("+")]
	elif "-"  in cdna_pos:
		[cdna_pos, offset] = [int(i) for i in cdna_pos.split("-")]
		offset = -offset
	if cdna_pos in abca4_cdna2gdna:
		print(f"{cdna_pos} not found in abca4_cdna2gdna")
	else:
		print(f"reconstructing value for {cdna_pos} in abca4_cdna2gdna")
	exit()


def parse_cdna(cdna_var):
	location = None
	mod_from = None
	mod_to = None
	[pos_from, pos_to] = [None, None]
	for mt in ["delins", "del", "ins", "dup"]:
		if not mt in cdna_var: continue
		mod_type = mt
		fields = cdna_var.split(mt)
		location = fields[0]
		if len(fields)>1: mod_to = fields[1]
		break
	if location:
		if "_" in location:
			[pos_from, pos_to] = location.split("_")
		else:
			[pos_from, pos_to] = [location, None]
	else: # this is substitution
		mod_type = "sub"
		pattern = re.match('([\d\-+]+)(\D)>(\D)', cdna_var)
		if not pattern:
			print(f"!! {cdna_var}")
			exit()
		pos_from = pattern.group(1)
		mod_from = pattern.group(2)
		mod_to   = pattern.group(3)
	print(f"\ncdna_var:{cdna_var}\n\tmod_type: {mod_type}\n\tpos_from: {pos_from}", end="")
	print(f"\n\tpos_to: {pos_to}\n\tmod_from: {mod_from}\n\tmod_to: {mod_to}")
	gdna_start = gdna_from_cdna(pos_from)
	gdna_end = gdna_from_cdna(pos_to) if pos_to else gdna_start
	return [gdna_start, gdna_end, mod_type, mod_from, mod_to]

def  store_variants(cursor, c, p):
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
		# print(f"storing {cdna_var} {protein_var}")
		if cdna_var!="np":
			[gdna_start, gdna_end, mod_type, nt_from, nt_to] = parse_cdna(cdna_var)
			qry = f"select * from variants where cdna='{cdna_var}'"
			ret = error_intolerant_search(cursor, qry)
			if ret:
				# check that we agree about the protein variant
				pass
			else:
				# qry  = "insert into variants (gdna_start, gdna_end, mod_type, nt_from, nt_to, protein) values "
				# qry += f"values ({gdna_start}, {gdna_end}, '{mod_type}', '{nt_from}'', '{nt_to}', '{protein}')"
				# print(qry)
				pass

		else: # see if we have the protein variant
			pass
			#print("we have an np")

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
		# store_or_update(cursor, "cases", fixed_fields={'publication_id': publication_id, 'patient_xref_id': patient_id},
		# 						update_fields ={'allele_id_1': allele_ids[0], 'allele_id_2': allele_ids[1],
		# 										'onset_age': onset, 'acuity_type':'BCVA', 'eye':'better', 'progression':progression_string})

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
