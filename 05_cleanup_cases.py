#!/usr/bin/python3

import os
import re

from utils.mysql import *
from utils.annotation import *
from sys import argv, stdout


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
perception_methods = ["finger", "hand", "fc", "hw", "light", "hm"]


def normalize(value_type, acuity):
	value_type = value_type.lower()
	for method in perception_methods:
		if method in value_type: return 0.0

	if "/" in acuity: # this is BCVA
		[nom,denom] = acuity.split("/")
		return float(nom)/float(denom)

	elif value_type.lower()=="logmar":
		try:
			# print(acuity, pow(10, -float(acuity)))
			return pow(10, -float(acuity))
		except:
			print (f"logMAR conversion failed for {value_type} {acuity}")
			exit()
	else:
		acuity = float(acuity)
		if acuity<0 or acuity>1.0:
			print(f"acuity out of range {value_type} {acuity}")
			exit()
	return 0.0


def progression_cleanup(progression, value_type):
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
	return ";".join(acuity_age)

#################################################################
def cdna_variant_list(cdna_allele_string):
	if not cdna_allele_string: return []
	# make sure we don't have the funny dash symbols for minus
	# - ord 45 ; − ord 8722 ; – ord 8211
	# the proper minus sign is 45
	cdna_allele_string = re.sub('[−–]', '-', cdna_allele_string)
	# no brackets
	cdna_allele_string = re.sub('[\[\]\(\)\s]', '', cdna_allele_string)\
		.replace("c.", "").replace("np", "").replace("NP", "").upper()\
		.replace("DUP","dup").replace("DEL","del").replace("INS","ins").replace("DELINS","delins")\
		.split(";")

	return list(filter(lambda x: x!=None and x!="", cdna_allele_string))

def protein_variant_list(protein_allele_string):
	# the usage of * and Ter is invonsistent -
	# one would expect * to be used in dingle letter rep, and Ter in three letter, but ut is not the case
	protein_allele = re.sub('[\[\]\(\)\s]', '', protein_allele_string).replace("p.", "").split(";")
	protein_allele = list(filter(lambda x: x!=None and x!="", protein_allele))
	return protein_allele

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

def aa2single(aa, protein_variant):
	if len(aa)==3:
		if aa.lower()=='del': return"del"
		aa_from = single_letter_code.get(aa.upper(), "none")
		if aa_from != "none": return aa_from
	elif len(aa)==1:
		if aa.upper() in single_letter_code.values(): return aa.upper()
	print(f"unrecognized aa: {protein_variant} ({aa}?)")
	exit()
	return aa


def protein_variant_sanity_check(protein_variant):

	if "splice" in protein_variant.lower(): return protein_variant
	if "deep" in protein_variant.lower(): return  protein_variant

	pattern1 = re.match('(\D+)(\d+)(\D+)', protein_variant)
	pattern2 = re.match('(\D+)(\d+)_(.+)', protein_variant)
	aa_from = None
	aa_to = None
	pos = -1
	fs_pos = None
	if pattern2:
		aa_from = pattern2.group(1)
		aa_from = aa2single(aa_from, protein_variant)
		pos     = pattern2.group(2)
		anything =  pattern2.group(3)
	elif pattern1:
		aa_from = pattern1.group(1)
		aa_from = aa2single(aa_from, protein_variant)
		pos     = pattern1.group(2)
		aa_to = pattern1.group(3)
		if "fs" in aa_to:
			print(aa_to)
			[aa_to, fs_pos] = aa_to.split("fs")
			fs_pos = fs_pos.replace("*", "Ter")
		aa_to = aa2single(aa_to, protein_variant)
	else:
		print(f"unrecognized protein variant: {protein_variant}")
		exit()

	protein_sequence =  get_protein()
	if int(pos)-1>=len(protein_sequence):
		print(f"protein length out of range for {protein_variant}")
		exit()

	orig_aa = protein_sequence[int(pos)-1]
	if aa_from != orig_aa:
		print(f"type mismatch for {protein_variant} the original AA type at {pos} is {orig_aa}, here is {aa_from}")
		exit()

	if pattern2:
		protein_variant = f"{three_letter_code[aa_from]}{pos}_{anything}"
	else:
		if fs_pos:
			protein_variant = f"{three_letter_code[aa_from]}{pos}{three_letter_code[aa_to]}fs{fs_pos}"
		elif aa_to in ['del']:
			protein_variant = f"{three_letter_code[aa_from]}{pos}{aa_to}"
		else:
			protein_variant = f"{three_letter_code[aa_from]}{pos}{three_letter_code[aa_to]}"

	return protein_variant


def allele_cleanup(cdna_allele_string, protein_allele_string,  verbose=False, outf=None):

	cdna_allele = cdna_variant_list(cdna_allele_string)
	protein_allele = protein_variant_list(protein_allele_string)
	lp = len(protein_allele)
	lc = len(cdna_allele)

	if lp==0 and lc==0: return [],[]

	if lp>0 and lc>0 and lp!=lc:
		print(f"differing number of variants in cDNA and protein description: {protein_allele} {cdna_allele}")
		exit()

	protein_allele_clean = []
	cdna_allele_clean = []
	if lp==0:
		for i in range(len(cdna_allele)):
			cdna_variant    = cdna_allele[i]
			cdna_variant_clean, protein_variant_inferred = mutation_effect(cdna_variant)
			if "splice" in protein_variant_inferred:
				if verbose: print(" %-20s  %-20s"%(cdna_variant, protein_variant_inferred), file=outf)
			else:
				if verbose: print(" %-20s  %-20s (inferred)"%(cdna_variant, protein_variant_inferred), file=outf)
			protein_allele_clean.append(protein_variant_inferred)
			cdna_allele_clean.append(cdna_variant_clean)

	elif lc==0:
		for i in range(len(protein_allele)):
			cdna_variant_clean = "np"
			protein_variant = protein_allele[i]
			protein_variant = protein_variant_sanity_check(protein_variant) # this will exit if not acceptable
			if verbose: print(" %-20s  %-20s"%(cdna_variant_clean, protein_variant), file=outf)
			protein_allele_clean.append(protein_variant)
			cdna_allele_clean.append(cdna_variant_clean)

	else:
		for i in range(len(cdna_allele)):
			cdna_variant    = cdna_allele[i]
			protein_variant = protein_allele[i]
			cdna_variant_clean, protein_variant_inferred = mutation_effect(cdna_variant)
			check = compare_protein_variants(protein_variant_inferred, protein_variant)
			if verbose: print(" %-20s  %-20s  %-20s    %s"%(cdna_variant, protein_variant, protein_variant_inferred, check), file=outf)
			protein_allele_clean.append(protein_variant_inferred)
			cdna_allele_clean.append(cdna_variant_clean)

	return cdna_allele_clean, protein_allele_clean


def parse_in(in_tsv, verbose=False, outf=None, faux=False):
	cases = {}
	inf = open(in_tsv)
	linect = 0
	reference = {}
	for line in inf:
		if verbose: print(line.strip())
		fields = [f.strip() for f in line.split('\t')]
		if 'pubmed' in fields[0].lower(): continue # this is the header
		linect +=1
		[pmid, ref, patient_id, cdna1, protein1, cdna2, protein2, hapl_tested, comment, value_type, onset] = fields[:11]
		reference[pmid] = ref.replace(" ","")
		try:
			onset = int(onset)
		except:
			onset = -1

		if faux:
			progression_string = ""
		else:
			progression = [f for f in fields[11:]]
			if len(progression)==0:
				if verbose: print("progression not given", pmid, ref, patient_id, file=outf)
				progression = None

			if progression and len(progression)%3!=0:
				print("length of progression not divisible by 3", pmid, ref, patient_id)
				exit()
			progression_string = progression_cleanup(progression, value_type) if progression else ""

		# if hapl_tested.lower() != 'yes':continue # let's stay away from those for now


		cdna_allele_1, protein_allele_1 = allele_cleanup(cdna1, protein1, verbose=verbose, outf=outf)
		if len(cdna_allele_1)==0 and len(protein_allele_1)==0:
			print("one allele missing in\n", line)
			exit()
		cdna_allele_2, protein_allele_2 = allele_cleanup(cdna2, protein2, verbose=verbose, outf=outf)
		if len(cdna_allele_2)==0 and len(protein_allele_2)==0:
			print("one allele missing in\n", line)
			exit()


		patient_key = f"{pmid}\t{patient_id}"
		if patient_key in cases:
			print("duplicate patient key:", patient_key)
			exit()
		# value type is BCVA, that is why we went through the "normalization" above
		cases[patient_key] = [";".join(cdna_allele_1), ";".join(protein_allele_1),
										";".join(cdna_allele_2), ";".join(protein_allele_2),
										'BCVA', onset, progression_string]

	inf.close()
	if verbose: print("lines read in", linect, file=outf)
	return cases, reference



#########################################
def main():
	# run this once with verbose True option to make sure you can live with all mismatches in the mutation interpretation
	# (or delete from input if you cannot)

	if len(argv)<2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()
	if not in_tsv[-4:]==".tsv":
		print(f"is {in_tsv} a tsv file?")
		exit()

	# logf = open("load_cases.log", "w")
	logf = None
	cases, reference = parse_in(in_tsv, verbose=False, outf=logf)

	out_tsv = in_tsv[:-4] + ".clean.tsv"
	with open(out_tsv, "w") as outf:
		for case, data in cases.items():
			[pubmed_id, patient_id] = case.split("\t")
			[c1, p1, c2, p2, value, onset, progression_string] = data
			print("\t".join([pubmed_id, reference[pubmed_id], patient_id, c1, p1, c2, p2,  "blah", "blah", value, str(onset), progression_string]),  file=outf)
	print("number of cases", len(cases), file=logf)

	# self consistency check
	# cases = parse_in(out_tsv, verbose=True, outf=logf, faux=True)




#########################################
if __name__ == '__main__':
	main()

'''
various experiments


	# last_nucleoutide_bedore_donor() # this really seems to be overwhelmingly G
	# first_nucleoutide_after_acceptor() $ this is more commonly G, but not as overwhlmingly as in the case of donor
	# exit()

	# codon = get_codons()[2029]
	# print(codon)
	# print(f"CGA  {Seq('CGA').translate()}   AGA {Seq('AGA').translate()}   TGA {Seq('TGA').translate()} ")
	# exit()

	#  5161_5162delAC and 5160_5161delCA have the same consequence on the protin level
	# print(get_cdna()[5159:5162])
	# codon = get_codons()
	# protein = get_protein()
	# for aa in range(1718,1726):
	# 	print(aa, aa*3, codon[aa-1], protein[aa-1], Seq(codon[aa-1]).translate())
	# print()
	
	# CAC
	# 1718 5154 GTG V V
	# 1719 5157 AGC S S
	# 1720 5160 CCC P P
	# 1721 5163 ACC T T
	# 1722 5166 ACC T T
	# 1723 5169 TAC Y Y
	# 1724 5172 TGG W W
	# 1725 5175 GTG V V
	# exit()



'''