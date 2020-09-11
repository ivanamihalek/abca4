#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from sys import argv, stdout
from utils.annotation import three_letter_code
from utils.abca4_gene import get_protein
from utils.utils import panic
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
		if len(fields)<1: continue
		if 'variant' in fields[0]: continue # this is header
		if len(fields)<9: continue # empty line
		variants.append(fields[:7]) # not interested in the last two fields
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

def store_publication(cursor, url, pmc, pubmed_id, ref):

	publication_id = None
	ret = error_intolerant_search(cursor, f"select id from publications where other_xref='{url}'")
	if not ret:
		qry = f"insert into publications (reference, other_xref) values ('{ref}', '{url}')"
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

def simplified(x):
	return float("%.1f"%(float(x.split("+/-")[0])/100))

#########################################
def process_input_tsv(cursor, tsv, reference):

	[url, pmc, pubmed_id, ref] =  reference
	variants = parse_in(tsv)
	protein_sequence =  get_protein()
	wt_values = []
	for variant in variants:
		[v, localization, solubilization, basal_atpase_activity,
		 n_ret_p_atpase_activity,  basal_cargo_binding, atp_cargo_binding] = variant
		if v=='WT':
			wt_values = [simplified(x) for x in variant[3:]]
			continue
		frm = v[0]
		to  = v[-1]
		pos = int(v[1:-1])
		protein_variant = f"{three_letter_code[frm]}{pos}{three_letter_code[to]}"
		# quick sanity check
		sanity(frm, pos, to, protein_sequence, protein_variant)
		variant_ids = store_variant(cursor, protein_variant)
		expression = simplified(solubilization)
		if expression>1: expression=1.0
		transport_values = [simplified(x) for x in variant[3:]]
		transport = 0
		for i in range(len(transport_values)):
			d = transport_values[i]/wt_values[i] if transport_values[i]<wt_values[i]  else 1.0
			transport += d*d
		transport = float("%.1f"%(transport/len(transport_values)))
		# print(frm, pos, to, expression, transport, url, pmc, pubmed_id, protein_variant, variant_ids)
		#
		# store or retrieve publication id
		publication_id = store_publication(cursor, url, pmc, pubmed_id, ref)
		for var_id in variant_ids:
			fixed_fields = {'variant_id':var_id}
			update_fields = {'expression_folding_membrane_incorporation':expression,
			                 'transport_efficiency':transport,
			                 'publication_id':publication_id}

			store_or_update(cursor, "parametrization_literature", fixed_fields=fixed_fields, update_fields=update_fields)


#########################################
def main():

	# tsvs = ["raw_data/molday_tm1.tsv", "raw_data/molday_tm2.tsv"]
	# [url, pmc, pubmed_id, ref] = ["https://www.biorxiv.org/content/10.1101/2020.08.28.272914v1.full",
	#                               None, None, "Garces2020BioRX"]

	# tsvs = ["raw_data/molday_random.tsv"]
	# [url, pmc, pubmed_id, ref] = ["https://iovs.arvojournals.org/article.aspx?articleid=2680974",
	#                               None, None, "Garces2018IOVS"]

	tsvs = ["raw_data/molday_N965S.tsv"]
	[url, pmc, pubmed_id, ref] = ["https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5886264",
	                              None, None, "Molday2018HMG"]


	db, cursor = abca4_connect()
	for tsv in tsvs:
		process_input_tsv(cursor, tsv, [url, pmc, pubmed_id, ref])
	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
