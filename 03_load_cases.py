#!/usr/bin/python3

import os
from utils.mysql import *
from sys import argv

# really bad
# finger counting, hand waving, light perception
perception_methods = ["finger", "hand", "fc", "hw", "light"]

def normalize(value_type, acuity):

	print("implement logmar recalc")
	exit()

	value_type = value_type.lower()
	for method in perception_methods:
		if method in value_type: return 0.0

	if "/" in acuity:
		[nom,denom] = acuity.split("/")
		return float(nom)/float(denom)
	else:
		return float(acuity)


def cdna_cleanup(cdna):
	print("implement cdna cleanup")
	exit()
	return ""


def cdna_protein(protein):
	print("implement cdna cleanup")
	exit()
	return ""


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

		if hapl_tested.lower() != 'yes':
			continue # let's stay away from those for now

		acuity_age = []
		for [age, od, os] in [progression[i:i+3] for i in range(0, len(progression), 3)]:
			print(pmid, ref, patient_id, age, od, os)
			if "/" in value_type:
				value_type_od, value_type_os = value_type.split("/")
			else:
				value_type_od, value_type_os = value_type, value_type
			od = -1 if od == '' else normalize(value_type_od, od)
			os = -1 if os == '' else normalize(value_type_os, os)
			if os==-1 and od==-1: continue
			acuity_age.append("%.1f %.2f"%(float(age),max(od,os)))
		try:
			onset = int(onset)
		except:
			onset = -1

		cdna1 = cdna_cleanup(cdna1)
		cdna2 = cdna_cleanup(cdna2)
		protein1 = protein_cleanup(protein1)
		cdna2 = cdna_cleanup(cdna2)

		cases[f"{pmid} {patient_id}"] = [cdna1, protein1, cdna2, protein2, value_type, onset, ";".join(acuity_age)]
	inf.close()
	return cases

#########################################
def main():

	if len(argv)<2:
		print(f"usage: {argv[0]} <input tsv>")
		exit()

	in_tsv = argv[1]
	if not os.path.exists(in_tsv):
		print(f"{in_tsv} not found")
		exit()

	cases = parse_in(in_tsv)
	for case, data in cases.items():
		print(case, data)

	db,cursor = abca4_connect()


	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

