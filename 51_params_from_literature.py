#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from sys import argv, stdout
from utils.annotation import three_letter_code
from utils.abca4_gene import get_protein
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein


def parse_in(in_tsv):
	variants = []
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()
	if not in_tsv[-4:] == ".tsv":
		print(f"is {in_tsv} a tsv file?")
		exit()
	inf = open(in_tsv)
	for line in inf:
		fields = line.strip().split('\t')
		if len(fields)<1:continue
		if 'frm' in fields[0]: continue # this is header
		if len(fields)<6:continue
		if len(fields)<8:
			fields += [None]*(8-len(fields))
		variants.append(fields)
	return variants


def store_variant(cursor, protein_var, verbose=False):

	qry = f"select  id from variants where protein='{protein_var}'"
	ret = error_intolerant_search(cursor, qry)
	if ret:
		# print(protein_var)
		# here we will make an assumption that the variant on the protein
		# level is causes by the same nucleotide change as the known one
		return [r[0] for r in ret]
	else:
		mod_type = "sub"
		qry = "insert into variants (mod_type, protein) "
		qry += f"values ('{mod_type}','{protein_var}')"
		if search_db(cursor, qry, verbose=verbose): return False
		return [hard_landing_search(cursor, "select max(id) from variants")[0][0]]

	return True

def sanity(frm, pos, to, protein_sequence, protein_variant):
	if int(pos)-1>=len(protein_sequence):
		print(f"protein length out of range for {protein_variant}")
		exit()
	protein_variant = f"{three_letter_code[frm]}{pos}{three_letter_code[to]}"
	orig_aa = protein_sequence[int(pos)-1]
	if frm != orig_aa:
		print(f"type mismatch for {protein_variant} the original AA type at {pos} is {orig_aa}, here is {frm}")
		exit()

def store_publication(cursor, url, pmc, pubmed_id):

	publication_id = None
	ret = error_intolerant_search(cursor, f"select id from publications where other_xref='{url}'")
	if not ret:
		qry = f"insert into publications (other_xref) values ('{url}')"
		if search_db(cursor, qry, verbose=True): exit()
		publication_id = hard_landing_search(cursor, "select max(id) from publications")[0][0]
		if pubmed_id:
			qry = f"update publications set pubmed={pubmed_id} where id={publication_id}"
			if search_db(cursor, qry, verbose=True): exit()
		if pmc:
			qry = f"update publications set pubmedcentral='{pmc}' where id={publication_id}"
			if search_db(cursor, qry, verbose=True): exit()

	elif len(ret)>1:
		panic(["multiple returns for", pubmed_id])
	else:
		publication_id = ret[0][0]
	return publication_id

#########################################
def main():

	if len(argv) < 2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	variants = parse_in(in_tsv)
	protein_sequence =  get_protein()

	db, cursor = abca4_connect()
	for variant in variants:
		[frm, pos, to, expression, transport, url, pmc, pubmed_id] = variant
		protein_variant = f"{three_letter_code[frm]}{pos}{three_letter_code[to]}"
		# quick sanity check
		sanity(frm, pos, to, protein_sequence, protein_variant)
		variant_ids = store_variant(cursor, protein_variant)
		print(frm, pos, to, expression, transport, url, pmc, pubmed_id, protein_variant, variant_ids)
		#
		# store or retrieve publication id
		publication_id = store_publication(cursor, url, pmc, pubmed_id)
		exit()
		#
		# #store case: allele_id_1, allele_id_2, publication_id, patient_id, onset, value, better, progression
		# store_or_update(cursor, "cases", fixed_fields=fixed_fields, update_fields=update_fields)
		#
		# # print(fixed_fields)
		# # print(update_fields)
		# # print()
		# #exit()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
