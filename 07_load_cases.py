#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from sys import argv, stdout
from utils.abca4_gene import *
from Bio.Seq  import Seq
from Bio.Alphabet import generic_dna, generic_protein


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


sorted_exon_bdries = sorted(abca4_cdna2gdna.keys())


###########################################################
def gdna_from_cdna(cdna_pos, verbose=False):
	# return the position, and the nucleotide, for sanity checking
	if cdna_pos is None:
		print("cdna cannot be none here (do something)")
		exit()
	offset = 0
	if "+" in cdna_pos:
		[cdna_pos, offset] = [int(i) for i in cdna_pos.split("+")]
	elif "-"  in cdna_pos:
		[cdna_pos, offset] = [int(i) for i in cdna_pos.split("-")]
		offset = -offset
	else:
		cdna_pos = int(cdna_pos)
	gdna_pos = None
	if cdna_pos in abca4_cdna2gdna:
		# with abca4 we are on the - strand
		offset = -offset
		gdna_pos = abca4_cdna2gdna[cdna_pos] + offset
		if verbose: print(f"{cdna_pos} is exon boundary;  gdna_pos = {gdna_pos}, offset = {offset}")
	else:
		if verbose: print(f"reconstructing value for {cdna_pos} in abca4_cdna2gdna")
		for i in range(1, len(sorted_exon_bdries)):
			p_prev = sorted_exon_bdries[i-1]
			p  = sorted_exon_bdries[i]
			if cdna_pos<p:
				# gdna_pos = abca4_cdna2gdna[p_prev] + cdna_pos - p_prev
				# we ar on the - strand
				gdna_pos = abca4_cdna2gdna[p_prev] - (cdna_pos - p_prev)
				if verbose: print(f"{p_prev} < {cdna_pos} < {p}  {abca4_cdna2gdna[p_prev]}  {abca4_cdna2gdna[p]}        {gdna_pos}")
				break
	if not gdna_pos:
		print("failed to reconstruct the gdna pos")
		exit()
	return gdna_pos


def get_nt_from_gdna(region_start, region_end):
	# this is a total hack that works only fo the current setup and chromosome 1
	# blastdbcmd, cdna_fasta, seq_region_name
	blastdbcmd = "/usr/bin/blastdbcmd"
	cdna_fasta = "/storage/databases/ensembl-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.fa"
	seq_region_name = "1"

	[frm, to] = [region_start, region_end] if  region_start <= region_end else [region_end, region_start]

	tmpfile = "tmp.fasta"
	if os.path.exists(tmpfile): os.remove(tmpfile)
	cmd  = f"{blastdbcmd} -db {cdna_fasta} -dbtype nucl -entry {seq_region_name} "
	cmd += f"-range {frm}-{to} -out {tmpfile} -outfmt %s"
	subprocess.call(["bash", "-c", cmd])
	if not os.path.exists(tmpfile):
		print(f"{tmpfile} not produced")
		exit()
	with open(tmpfile) as inf:
		inseq = inf.read().replace("\n", "")
	return inseq


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
	gdna_start = gdna_from_cdna(pos_from)
	gdna_end = gdna_from_cdna(pos_to) if pos_to else gdna_start

	# sanity
	gdna_nt = Seq(get_nt_from_gdna(gdna_start, gdna_end)).reverse_complement()

	if mod_type not in ["del","ins", "delins"] and mod_from!=gdna_nt:
		print(f"\ncdna_var:{cdna_var}\n\tmod_type: {mod_type}\n\tpos_from: {pos_from}", end="")
		print(f"\n\tpos_to: {pos_to}\n\tmod_from: {mod_from}  gdna_nt: {gdna_nt} \n\tmod_to: {mod_to}", end="")
		print(f"\n\tgdna_start: {gdna_start}\n\tgdna_end: {gdna_end}")
		exit()
	if mod_type=="del":
		mod_from = gdna_nt
		mod_to = "-"
	elif mod_type=="ins":
		mod_from = "-"
	elif  mod_type=="delins":
		mod_from = gdna_nt

	return [gdna_start, gdna_end, mod_type, mod_from, mod_to]


###############################################
def store_variants(cursor, c, p, verbose=True):
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

	var_ids = []
	for i in range(len(cdna_vars)):
		cdna_var = cdna_vars[i].replace(" ","")
		protein_var = protein_vars[i].replace(" ","")
		if cdna_var=="np":
			if not protein_var or protein_var=="np":
				print("neither cdna nor protein provided", cdna_vars, protein_vars)
				exit()
			qry = f"select  id from variants where protein='{protein_var}'"
			ret = error_intolerant_search(cursor, qry)
			if ret:
				# print(protein_var)
				# here we will make an assumption that the variant on the protein
				# level is causes by the same nucleotide change as the known one
				if len(ret)>1:
					print("    found", len(ret), "not sure what to assign to")
					print(qry)
					print(ret)
					exit()
				var_ids.append(ret[0][0])
			else:
				mod_type = "sub"
				if "fs" in protein_var:
					# on the protein level we cannot know how the frameshift arose
					# perhaps  it would not be that hard to guess?
					mod_type = "unk"
				else:
					for mt in ["delins", "del", "ins", "dup"]:
						if mt in protein_var:
							mod_type = mt
							break
				qry  = "insert into variants (mod_type, protein) "
				qry += f"values ('{mod_type}','{protein_var}')"
				if search_db(cursor, qry, verbose=verbose): exit()
				var_ids.append(hard_landing_search(cursor, "select max(id) from variants")[0][0])

		else:
			[gdna_start, gdna_end, mod_type, nt_from, nt_to] = parse_cdna(cdna_var)
			qry = f"select  id, gdna_start, gdna_end, protein from variants where cdna='{cdna_var}' or protein='{protein_var}'"
			ret = error_intolerant_search(cursor, qry)
			if ret:
				# this is an existing entry - just check that we agree on the protein variant
				if len(ret)>1:
					print(f"more than one entry present for {cdna_var}", ret)
					exit()
				if verbose: print(f"storing {cdna_var} found:", ret[0])
				[var_id, gdna_start_in_db, gdna_end_in_db, protein_in_db] = ret[0]
				# TODO - if the stored gdna is null, than we already bumped into this
				# variant without cdna provided - use this entry to fill
				if gdna_start_in_db != gdna_start or gdna_end_in_db != gdna_end or protein_in_db != protein_var:
					print(f"mismatch for {cdna_var}: ")
					print(gdna_start_in_db, gdna_end_in_db, protein_in_db)
					print(gdna_start, gdna_end, protein_var)
					exit()
				var_ids.append(var_id)
			else:
				if verbose: print(f"storing {cdna_var} {protein_var}")
				if len(nt_from)>10: nt_from = nt_from[:5]+"..."+ nt_from[-5:]
				if len(nt_to)>10: nt_to = nt_to[:5]+"..."+ nt_to[-5:]
				qry  = "insert into variants (gdna_start, gdna_end, mod_type, nt_from, nt_to, cdna, protein) "
				qry += f"values ({gdna_start}, {gdna_end}, '{mod_type}', '{nt_from}', '{nt_to}', '{cdna_var}', '{protein_var}')"

				if search_db(cursor, qry, verbose=verbose): exit()
				var_ids.append(hard_landing_search(cursor, "select max(id) from variants")[0][0])
				#print()

	if len(var_ids)!= len(cdna_vars):
		print(f"missing variant id for {cdna_vars} {protein_vars} (?)")
		exit()
	return var_ids

#########################################
def main():


	# parse_cdna('972_973delinsAT')
	# # #print(get_cdna()[3810])
	# # #parse_cdna('4510_4535del')
	# # #print(get_cdna()[4509:4535])
	# exit()

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
			variant_ids = store_variants(cursor, cdna_variants, protein_variants)
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

