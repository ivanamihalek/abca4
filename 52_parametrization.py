#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.structure import *
from utils.utils import  *

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
def annotate_by_cons(cons_data, protein, protein_domain):
	[expression, transport, notes] = [1.0, 1.0, ""]
	[type_in_specs, substitutions, cons] = cons_data
	[aa_from_long, pos, aa_to_long] = parse_protein(protein)
	aa_from = single_letter_code[aa_from_long]
	aa_to = single_letter_code[aa_to_long]
	if aa_from!=type_in_specs[pos]:
		print(f"type mismatch for pos {pos}: {aa_from} vs {type_in_specs[pos]} in specs")
		exit()
	if aa_to in substitutions[pos]:
		notes = "mild: substitution seen"
	else:
		if protein_domain in ["NBD1", "NBD2"]:
			transport = 0.70
			notes = "transport: NBD, novel subs"
		elif protein_domain in ["TMD1", "TMD2"]:
			transport = 0.50
			notes = "transport: TMD, novel subs"
		elif protein_domain in ["ECD1", "ECD2"]:
			transport = 0.90
			notes = "transport: ECD, novel subs"
		else:
			transport = 0.90
			notes = f"transport: {protein_domain}, novel subs"

	return [expression, transport, notes]


#########################################
def main():

	cons_data = read_specs("conservation/vertebrates_only/specs_out.score")

	db, cursor = abca4_connect()
	# find column names
	qry  = "select column_name from information_schema.columns "
	qry += "where table_schema='abca4' and  table_name='variants' order by ordinal_position"
	header = [row[0] for row in hard_landing_search(cursor, qry)]

	total = 0
	annotated = 0
	for fields in hard_landing_search(cursor, "select * from variants"):
		total += 1
		annotated += 1
		named_field = dict(zip(header, fields))
		protein = named_field["protein"]
		cdna = named_field["cdna"]
		pd = named_field["protein_domain"]
		# initial guess: both factors unaffected
		[expression, transport, notes] = [1.0, 1.0, ""]
		if "fs" in protein:
			expression = 0.0
			notes = "null: fs"
		elif is_splice(cdna, protein):
			expression = 0.0
			notes = "null: near splice"
		elif is_ter(cdna, protein):
			expression = 0.0
			notes = "null: early stop"
		elif is_del(cdna, protein):
			expression = 0.0
			notes = "null: del"
		elif is_exotic(cdna, protein):
			expression = 0.0
			notes = "null: exotic"
		elif is_misfolder(cdna, protein):
			expression = 0.0
			notes = "null: misfolder"
			#print(protein, "misfolder ===>", named_field["conserved_in_verts_insects"])
		elif is_distant_splice(cdna, protein):
			expression = 0.5
			notes = "expression: distant splice"

		elif is_synonymous(protein):
			notes = "suspicious: synonymous"

		elif named_field["conserved_in_verts_insects"]==1:
			expression = 0.0
			notes = "null: cons in para"

		elif named_field["conserved_in_ortho_verts"]==1:
			if pd in ["NBD1", "NBD2"]:
				transport = 0.50
				notes = "transport: NBD, cons in ortho"
			elif pd in ["TMD1", "TMD2"]:
				transport = 0.25
				notes = "transport: TMD cons in ortho"
			elif pd in ["ECD1", "ECD2"]:
				transport = 0.75
				notes = "transport: ECD, cons in ortho"
			else:
				expression = 0.80
				transport = 0.80
				notes = f"expr and transp: {pd}, cons in ortho"


		elif "ins" in protein or "del" in  protein or "dup" in protein:
			if pd=="linker":
				expression = 0.5
				transport = 0.5
				notes = "null: medium structural mod"
			else:
				expression = 0.0
				notes = "null: gross structural mod"

		elif "deep" in protein :
			notes = "suspicious: putative deep intronic"

		elif is_missense(protein):
			[expression, transport, notes] = annotate_by_cons(cons_data, protein, pd)
		else:
			annotated -= 1
			print(cdna, protein)

		print("\t", expression, transport, notes)
		store_or_update(cursor, 'parametrization', fixed_fields={"variant_id":named_field["id"]},
		                update_fields={"expression_folding_membrane_incorporation":expression,
		                               "transport_efficiency":transport,
		                               "notes": notes})

	print(total, annotated)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()

