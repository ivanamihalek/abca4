#!/usr/bin/python3 -u


from utils.annotation import  single_letter_code
from utils.abca4_mysql import  *
import re
import matplotlib.pyplot as plt

def is_ter(cdna, prot):
	if "ter"in prot.lower(): return True
	return False

def is_del(cdna, prot):
	if "del"in prot.lower(): return True
	return False

# this is not very reliable:
# we have a case wiht two "misfolder" alleles that is over 40 at onset
def is_misfolder(cdna, protein_vars):
	if "_" in protein_vars: return False
	for protein in protein_vars.split(";"):
		pattern = re.findall('(\D{3})(\d+)(\D{3})', protein.strip().replace("p.",""))
		if not pattern: return False
		for match in pattern:
			aa_from = match[0].upper()
			pos = match[1]
			aa_to = match[2].upper()
			print(f"{protein}  {aa_from}  {aa_to} ")
			print(f"{single_letter_code[aa_from]}{pos}{single_letter_code[aa_to]}")
			if f"{single_letter_code[aa_from]}{pos}{single_letter_code[aa_to]}" in ["A1357T", "A1794P", "L2027F", "R2077W"]:
				return True
	return False


def is_exotic(cdna, protein):
	if "5461-10" in cdna: return True
	return False

def is_splice(cdna, protein):
	pattern = re.findall('(\d+)[\-+](\d+)(\D)', cdna)
	if pattern:
		min_match = 50
		for match in pattern:
			if int(match[1])<min_match: min_match=int(match[1])
		if min_match<3: return True

	# last nucleotide before splice
	pattern = re.findall('(\D+)(\d+)(\D+)splice', protein)
	if pattern: return True

	return False

# if "del"in prot.lower(): return True
# if "5461-10T>C" in cdna: return True # we have exp proof of this one
def is_null(cdna, prot, filter):
	if len(filter)==0: return True
	if "ter" in filter and is_ter(cdna, prot): return True
	if "del" in filter and is_del(cdna, prot): return True
	if "splice" in filter and is_splice(cdna, prot): return True
	if "exotic" in filter and is_exotic(cdna, prot): return True
	# if "misfolder" in filter and is_misfolder(cdna, prot): return True
	return False

def characterize(cdna, prot, homozygs, regions, cons):
	if is_null(cdna, prot, "ter|splice|del|exotic"): return "null"

	homs = homozygs.split(";")
	if len(homs)==1 and int(homs[0].strip())>1: return "homozg"

	regs =  regions.split(";")
	if len(regs)==1:
		#if cons=="n": return None
		for domain in ["nbd", "tm", "ecd"]:
			if domain in regs[0].lower(): return domain
	return None

#########################################
def main():


	db, cursor = abca4_connect()
	number_of_bins = 75
	bin_size = 1
	bins = []
	for b in range(number_of_bins):
		bins.append(b*bin_size+0.1)

	allele_type = ["null", "nbd", "tm", "ecd", "homozg"]
	ages = {}
	qry = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		[cdna1, prot1, freqs1, homozygs1, region1, cons1] = variants_from_allele(cursor, allele_id_1)
		[cdna2, prot2, freqs2, homozygs2, region2, cons2] = variants_from_allele(cursor, allele_id_2)
		if onset_age<0: continue # the onse age not given

		type1 = characterize(cdna1, prot1, homozygs1, region1, cons1)
		type2 = characterize(cdna2, prot2, homozygs2, region2, cons2)
		if not type1 or not type2: continue
		key = " ".join(sorted([type1, type2]))
		if not key  in ages: ages[key] = []
		ages[key].append(onset_age)

	db.close()

	tot = 0
	for types, onset_ages in ages.items():
		[type1, type2] = types.split()
		ct = len(onset_ages)
		tot += ct
		print(f"{type1}  {type2}   {ct}")
	print(tot)

	rows = 3
	cols = 5
	fig, axs = plt.subplots(rows, cols)
	fig.text(0.5, 0.04, 'onset age', ha='center')
	fig.text(0.04, 0.5, 'number of patients', va='center', rotation='vertical')
	all_type_pairs = list(ages.keys())
	for n in range(min(rows*cols, len(all_type_pairs))):
		name = all_type_pairs[n]
		axs.flat[n].hist(ages[name], bins)
		axs.flat[n].legend([f"{name} ({len(ages[name])})"])

	#axs.flat[rows * cols-1].hist(age_set["all"], bins, cumulative=True, histtype='step')

	plt.show()



#########################################
if __name__ == '__main__':
	main()
