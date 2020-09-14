#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.exp_characterization import *
from utils.utils import  *
from utils.rule_based_parametrization import rule_based_params

#########################################
def read_specs(infile):
	if not os.path.exists(infile):
		print(infile, "not found")
		exit()
	cons = {}
	type_in_specs = {}
	substitutions = {}
	inf = open(infile, "r")
	header = None
	for line in inf:
		if not header:
			# the first field is the comment sign
			header = line.strip().split()[1:]
		else:
			named_field = dict(zip(header,line.strip().split()))
			pos = int(named_field["pos_in_human"])
			type_in_specs[pos] = named_field["human"]
			substitutions[pos] = named_field["substitutions"]
			cons[pos] = named_field["rvet"]
	inf.close()

	return [type_in_specs, substitutions, cons]


#########################################
def annotate_variable_position(cons_data, protein, protein_domain):
	[expression, transport, notes] = [1.0, 1.0, ""]
	[type_in_specs, substitutions, cons] = cons_data
	[aa_from_long, pos, aa_to_long] = parse_protein(protein)
	aa_from = single_letter_code[aa_from_long]
	aa_to = single_letter_code[aa_to_long]
	if aa_from!=type_in_specs[pos]:
		print(f"type mismatch for pos {pos}: {aa_from} vs {type_in_specs[pos]} in specs")
		exit()
	if aa_to in substitutions[pos]:
		notes = "none: substitution seen"
	else:
		if protein_domain in ["NBD1", "NBD2"]:
			transport = 0.75
			notes = "mild transport: NBD, novel subs"
		elif protein_domain in ["TMD1", "TMD2"]:
			transport = 0.50
			notes = "strong transport: TMD, novel subs"
		elif protein_domain in ["ECD1", "ECD2"]:
			transport = 0.25
			notes = "strong transport: ECD, novel subs"
		elif protein_domain in ["R"]:
			transport = 0.50
			notes = "sever transport: R, novel subs"
		else:
			transport = 0.75
			notes = f"mild transport: {protein_domain}, novel subs"

	return [expression, transport, notes]

def find_experimentally_characterized_protein(cursor):
	experiment = {}
	qry   = "select  v.protein, p.expression_folding_membrane_incorporation, p.transport_efficiency "
	qry  += "from parametrization_literature p left join variants v on p.variant_id=v.id"
	ret = error_intolerant_search(cursor, qry)
	if ret:
		for line in ret:
			experiment[line[0]] = line[1:]
	return experiment

#########################################
def main():

	cons_data = read_specs("../conservation/abca4/vertebrates_only/specs_out.score")

	db, cursor = abca4_connect()
	# find column names
	qry  = "select column_name from information_schema.columns "
	qry += "where table_schema='abca4' and  table_name='variants' order by ordinal_position"
	header = [row[0] for row in hard_landing_search(cursor, qry)]

	caracterized_protein = find_experimentally_characterized_protein(cursor)

	total = 0
	annotated = 0
	for fields in hard_landing_search(cursor, "select * from variants"):
		total += 1
		annotated += 1
		named_field = dict(zip(header, fields))
		protein = named_field["protein"]
		cdna = named_field["cdna"]
		# initial guess: both factors unaffected
		[expression, transport, notes] = [1.0, 1.0, ""]

		############################################
		# first let's get the handful of values for which we have exp support  out of the way
		if cdna in characterized_cdna:
			[expression, transport, notes] = characterized_cdna_values(cdna)
		elif protein in caracterized_protein:
			[expression, transport] = caracterized_protein[protein]
		else:
			[expression, transport, notes] = rule_based_params(named_field, cons_data)
		############################################

		print("\t", expression, transport, notes)
		store_or_update(cursor, 'parametrization', fixed_fields={"variant_id":named_field["id"]},
		                update_fields={"expression_folding_membrane_incorporation":expression,
		                               "transport_efficiency":transport, "notes": notes})

	print(total, annotated)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()

