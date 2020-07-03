#!/usr/bin/python3

import os
import re

from utils.mysql import *
from utils.annotation import *
from sys import argv


#########################################
# some exploration: the sequence before
# the donor splice site convincingly end in G in the majority of cases

def last_nucleoutide_bedore_donor():
	cdna = get_cdna()
	count={}
	for nt in list("ACTG"): count[nt] = 0
	for pos in abca4_donor_splice.keys():
		if pos<1:continue
		print(pos, cdna[pos-1])
		count[cdna[pos-1]] += 1
	print(count)


def first_nucleoutide_after_acceptor():
	cdna = get_cdna()
	count={}
	for nt in list("ACTG"): count[nt] = 0
	for pos in abca4_acceptor_splice.keys():
		if pos>=len(cdna):continue
		print(pos, cdna[pos-1])
		count[cdna[pos-1]] += 1
	print(count)
#########################################


# really bad
# finger counting, hand waving, light perception
perception_methods = ["finger", "hand", "fc", "hw", "light"]


def normalize(value_type, acuity):
	value_type = value_type.lower()
	for method in perception_methods:
		if method in value_type: return 0.0

	if "/" in acuity:
		[nom,denom] = acuity.split("/")
		return float(nom)/float(denom)

	elif value_type.lower()=="logmar":
		try:
			# print(acuity, pow(10, -float(acuity)))
			return pow(10, -float(acuity))
		except:
			return 0.0
	return 0.0


def cdna_cleanup(cdna_variant):
	if not cdna_variant: return []
	cdna_variant = re.sub('[\[\]\(\)\s]', '', cdna_variant)\
		.replace("c.", "").upper()\
		.replace("DUP","dup").replace("DEL","del").replace("INS","ins").replace("DELINS","delins")\
		.split(";")
	# sometimes people are confused about the splice numbering, in particular
	# if the gene is encoded on - strand
	# see if we can fix that
	variant = [mutation_effect(cv) for cv in cdna_variant]
	return list(filter(lambda x: x!=None and x!="", variant))

def aa_equal(a1, a2):
	if a1==a2: return True
	if len(a1)==1 and len(a2)==3:
		if three_letter_code.get(a1, "none").upper()==a2.upper(): return True
	if len(a2)==1 and len(a1)==3:
		if three_letter_code.get(a2, "none").upper()==a1.upper(): return True
	return False

def compare_protein_variants(p1,p2):

	if p1==p2: return "check ok"
	# the fs without specifying where the Ter is downstream is considered and acceptable shorthand
	if "fs" in p1 and "fs" in p2:
		if p1.split("fs")[0] == p2.split("fs")[0]: return "check ok"
	if "splice" in p1 and "splice" in p2: return  "check ok"

	# do we have short/long aa symbol difference?
	pattern_1 =  re.match('(\D+)(\d+)(\D+)', p1)
	pattern_2 =  re.match('(\D+)(\d+)(\D+)', p2)
	if pattern_1 and pattern_2:
		aa_from_1 = pattern_1.group(1)
		aa_from_2 = pattern_2.group(1)
		pos_1 = pattern_1.group(2)
		pos_2 = pattern_2.group(2)
		aa_to_1 = pattern_1.group(3)
		aa_to_2 = pattern_2.group(3)
		if pos_1==pos_2 and aa_equal(aa_from_1, aa_from_2) and aa_equal(aa_to_1, aa_to_2): return "check ok"


	return "<===== mismatch"


def protein_variant_sanity_check(protein_variant):
	pattern = re.match('(\D+)(\d+)(\D+)', protein_variant)
	if not pattern: return "no pattern" # nothing much we can do here
	aa_from = pattern.group(1)
	if len(aa_from)==3: aa_from = single_letter_code.get(aa_from.upper(), "none")
	pos     = pattern.group(2)
	aa_to   = pattern.group(3)
	if aa_to.lower() in ['dup']: return "nostandard aa" # nothing much we can do here
	protein_sequence =  get_protein()
	if int(pos)-1>=len(protein_sequence):
		print(f"protein length out of range for {protein_variant}")
		exit()

	orig_aa = protein_sequence[int(pos)-1]
	if aa_from != orig_aa:
		print(f"type mismatch for {protein_variant} the original AA type at {pos} is {orig_aa}, here is {aa_from}")
		exit()
	return "check ok"

def protein_cleanup(protein_allele_string, cdna_allele, verbose=False):

	# the usage of * and Ter is invonsistent -
	# one would expect * to be used in dingle letter rep, and Ter in three letter, but ut is not the case
	protein_allele = re.sub('[\[\]\(\)\s]', '', protein_allele_string).replace("p.", "").split(";")
	protein_allele = list(filter(lambda x: x!=None and x!="", protein_allele))

	lp = len(protein_allele)
	lc = len(cdna_allele)

	if lp==0 and lc==0: return []

	if lp>0 and lc>0 and lp!=lc:
		print(f"differing number of variants in cDNA and protein description: {protein_allele} {cdna_allele}")
		exit()

	protein_allele_clean = []
	if lp==0:
		for i in range(len(cdna_allele)):
			cdna_variant    = cdna_allele[i]
			protein_variant = mutation_effect(cdna_variant)
			if "splice" in protein_variant:
				if verbose: print(" %-20s  %-20s"%(cdna_variant, protein_variant))
				protein_variant = "splice" # use as the place holder in case the allele has multiple variants
			else:
				if verbose: print(" %-20s  %-20s (inferred)"%(cdna_variant, protein_variant))
			protein_allele_clean.append(protein_variant)

	elif lc==0:
		for i in range(len(protein_allele)):
			cdna_variant = "(not given)"
			protein_variant = protein_allele[i]
			protein_variant_sanity_check(protein_variant)
			if verbose: print(" %-20s  %-20s"%(cdna_variant, protein_variant))
			protein_allele_clean.append(protein_variant)

	else:
		for i in range(len(cdna_allele)):
			cdna_variant    = cdna_allele[i]
			protein_variant = protein_allele[i]
			protein_variant_inferred = mutation_effect(cdna_variant)
			check = compare_protein_variants(protein_variant_inferred, protein_variant)
			if verbose: print(" %-20s  %-20s  %-20s    %s"%(cdna_variant, protein_variant, protein_variant_inferred, check))
			protein_allele_clean.append(protein_variant_inferred)
	return protein_allele_clean


def parse_in(in_tsv, verbose=False):
	cases = {}
	inf = open(in_tsv)
	linect = 0
	for line in inf:
		# print(line.strip())
		fields = [f.strip() for f in line.split('\t')]
		if 'pubmed' in fields[0].lower(): continue # this is the header
		linect +=1
		[pmid, ref, patient_id, cdna1, protein1, cdna2, protein2, hapl_tested, comment, value_type, onset] = fields[:11]
		progression = [f for f in fields[11:]]
		if len(progression)==0:
			if verbose: print("progression not given", pmid, ref, patient_id)
			progression = None

		if progression and len(progression)%3!=0:
			print("length of progression not divisible by 3", pmid, ref, patient_id)
			exit()

		# if hapl_tested.lower() != 'yes':
		# 	continue # let's stay away from those for now
		progression_string = ""
		if progression:
			acuity_age = []
			for [age, od, os] in [progression[i:i+3] for i in range(0, len(progression), 3)]:
				if "/" in value_type:
					value_type_od, value_type_os = value_type.split("/")
				else:
					value_type_od, value_type_os = value_type, value_type
				od = od.replace(" ","")
				os = os.replace(" ","")
				od = -1 if od == '' else normalize(value_type_od, od)
				os = -1 if os == '' else normalize(value_type_os, os)
				if os==-1 and od==-1: continue

				if len(age.replace(' ',''))==0:
					acuity_age.append("")
				else:
					acuity_age.append("%.1f:%.2f"%(float(age),max(od,os)))
			try:
				onset = int(onset)
			except:
				onset = -1
			progression_string = ";".join(acuity_age)

		cdna_allele_1 = cdna_cleanup(cdna1)
		cdna_allele_2 = cdna_cleanup(cdna2)
		protein_allele_1 = protein_cleanup(protein1,cdna_allele_1, verbose)
		protein_allele_2 = protein_cleanup(protein2,cdna_allele_2, verbose)

		patient_key = f"{pmid} {patient_id}"
		if patient_key in cases:
			print("duplicate patient key:", patient_key)
			exit()
		# value type is BCVA, that is why we went through the "normalization" above
		cases[patient_key] = [";".join(cdna_allele_1), ";".join(protein_allele_1),
										";".join(cdna_allele_2), ";".join(protein_allele_2),
										'BCVA', onset, progression_string]
	inf.close()
	print("line ct", linect)
	return cases



#########################################
def main():
	# run this once with verbose True option to make sure you can live with all mismatches in the mutation interpretation
	# (or delete from input if you cannot)

	# last_nucleoutide_bedore_donor() # this really seems to be overwhelmingly G
	# first_nucleoutide_after_acceptor() $ this is more commonly G, but not as overwhlmingly as in the case of donor
	# exit()

	# codon = get_codons()[2029]
	# print(codon)
	# print(f"CGA  {Seq('CGA').translate()}   AGA {Seq('AGA').translate()}   TGA {Seq('TGA').translate()} ")
	# exit()

	if len(argv)<2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()

	# todo: get verbose to log file
	cases = parse_in(in_tsv, verbose=False)
	print("number of cases", len(cases))
	# for case, data in cases.items():
	# 	print(case, data)


	db,cursor = abca4_connect()


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

