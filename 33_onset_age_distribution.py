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
	if "misfolder" in filter and is_misfolder(cdna, prot): return True
	return False


#########################################
def main():


	db, cursor = abca4_connect()
	number_of_bins = 75
	bin_size = 1
	bins = []
	for b in range(number_of_bins):
		bins.append(b*bin_size+0.1)
	filters = ["ter", "splice", "del", "exotic",  "ter|splice|del|exotic"]
	age_set  = {"all":[]}
	for name in filters: age_set[name] = []
	qry = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		cdna1, prot1, freqs, homs = variants_from_allele(cursor, allele_id_1)
		cdna2, prot2, freqs, homs = variants_from_allele(cursor, allele_id_2)
		if onset_age<0: continue # the onse age not given

		age_set["all"].append(onset_age)
		for fltr in filters:
			if is_null(cdna1, prot1, fltr) and is_null(cdna2, prot2, fltr):
				age_set[fltr].append(onset_age)

	cursor.close()
	db.close()

	rows = 3
	cols = 3
	fig, axs = plt.subplots(rows, cols)
	fig.text(0.5, 0.04, 'onset age', ha='center')
	fig.text(0.04, 0.5, 'number of patients', va='center', rotation='vertical')
	age_set_names = ["all"] + filters
	for n in range(len(age_set_names)):
		name = age_set_names[n]
		axs.flat[n].hist(age_set[name], bins)
		axs.flat[n].legend([f"{name} ({len(age_set[name])})"])

	axs.flat[rows * cols-1].hist(age_set["all"], bins, cumulative=True, histtype='step')

	plt.show()


#########################################
if __name__ == '__main__':
	main()
