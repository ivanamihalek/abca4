#!/usr/bin/python3

import os
import re

from utils.mysql import *
from utils.annotation import *
from sys import argv



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
	cdna_variant = re.sub('[\[\]\(\)\s]', '', cdna_variant).replace("c.", "").split(";")
	return list(filter(lambda x: x!=None and x!="", cdna_variant))

def compare_protein_variants(p1,p2):
	if p1==p2: return "check ok"
	# the fs without specifying where the Ter is doenstream is considered and acceptable shorthand
	if "fs" in p1 and "fs" in p2:
		if p1.split("fs")[0] == p2.split("fs")[0]: return "check ok"
	# do we have short/long aa symbol difference?
	if len(p1)-len(p2)==4:
		if p1 == (three_letter_code.get(p2[0],"") + p2[1:-1] + three_letter_code.get(p2[-1],"")): return "check ok"
	if len(p2)-len(p1)==4:
		if p2 == (three_letter_code.get(p1[0],"") + p1[1:-1] + three_letter_code.get(p1[-1],"")): return "check ok"


	return "<===== mismatch"


def protein_cleanup(protein_allele_string, cdna_allele):

	protein_allele = re.sub('[\[\]\(\)\s]', '', protein_allele_string).replace("p.", "").replace("*", "Ter").split(";")
	protein_allele = list(filter(lambda x: x!=None and x!="", protein_allele))

	lp = len(protein_allele)
	lc = len(cdna_allele)

	if lp==0 and lc==0: return []

	if lp>0 and lc>0 and lp!=lc:
		print(f"differing number of variants in cDNA and protein description: {protein_allele} {cdna_allele}")
		exit()

	if lp==0:
		for i in range(len(cdna_allele)):
			cdna_variant    = cdna_allele[i]
			protein_variant = mutation_effect(cdna_variant)
			if "splice" in protein_variant:
				print(" %-20s  %-20s"%(cdna_variant, protein_variant))
			else:
				print(" %-20s  %-20s (inferred)"%(cdna_variant, protein_variant))
	elif lc==0:
		for i in range(len(protein_allele)):
			cdna_variant = "(not given)"
			protein_variant = protein_allele[i]
			print(" %-20s  %-20s"%(cdna_variant, protein_variant))
	else:
		# TODO: check whether the provided and the inferred protein variant match
		for i in range(len(cdna_allele)):
			cdna_variant    = cdna_allele[i]
			protein_variant = protein_allele[i]
			protein_variant_inferred = mutation_effect(cdna_variant)
			check = compare_protein_variants(protein_variant_inferred, protein_variant)
			print(" %-20s  %-20s  %-20s    %s"%(cdna_variant, protein_variant, protein_variant_inferred, check))

	return []


def parse_in(in_tsv):
	cases = {}
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

		# if hapl_tested.lower() != 'yes':
		# 	continue # let's stay away from those for now

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

		cdna_allele_1 = cdna_cleanup(cdna1)
		cdna_allele_2 = cdna_cleanup(cdna2)
		protein_allele_1 = protein_cleanup(protein1,cdna_allele_1)
		protein_allele_2 = protein_cleanup(protein2,cdna_allele_2)

		cases[f"{pmid} {patient_id}"] = [";".join(cdna_allele_1), ";".join(protein_allele_1),
										";".join(cdna_allele_2), ";".join(protein_allele_2),
										value_type, onset, ";".join(acuity_age)]
	inf.close()

	return cases

#########################################
def main():

	# print(mutation_effect("5044_5058del15"))
	# exit()

	if len(argv)<2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()

	cases = parse_in(in_tsv)
	print("number of cases", len(cases))
	# for case, data in cases.items():
	# 	print(case, data)

	db,cursor = abca4_connect()


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

