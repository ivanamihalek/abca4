#!/usr/bin/python3 -u


from utils.abca4_mysql import  *
import re
import matplotlib.pyplot as plt


def is_null(cdna, prot):
	if "ter"in prot.lower(): return True
	if "del"in prot.lower(): return True
	# if "5461-10T>C" in cdna: return True # we have exp proof of this one
	# # match matches the beginning of the string only
	# # pattern = re.amtch('(\d+)[\-+](\d+)(\D)>(\D)', cdna)
	# # match function has groups, e.g pattern.group(2)
	# # the offset is 1, bacuse the group(0) is the overall match
	# # findall returns a list, and has 0 offset
	#
	# # mutations in the first two positions in the splice site
	# # also, insert or deletion that starts in the rist ro second position in the scplice site
	pattern = re.findall('(\d+)[\-+](\d+)(\D)', cdna)
	if pattern:
		min_match = 50
		for match in pattern:
			if int(match[1])<min_match: min_match=int(match[1])
		if min_match<3:  return True
		return False
	# the original splice pos not reported, or the last nulcotide befor the splice
	if "splice"in prot.lower(): return True
	return False


#########################################
def main():


	db, cursor = abca4_connect()
	number_of_bins = 64
	bin_size = 1
	bins = []
	for b in range(number_of_bins):
		bins.append(b*bin_size+0.1)
	age_set  = {}
	age_set_names = ["all", "both null", "one null", "none null"]
	for name in age_set_names: age_set[name] = []
	qry = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		cdna1, prot1 = variants_from_allele(cursor, allele_id_1)
		cdna2, prot2 = variants_from_allele(cursor, allele_id_2)
		if onset_age<0: continue
		age_set["all"].append(onset_age)
		if is_null(cdna1, prot1) and  is_null(cdna2, prot2):
			age_set["both null"].append(onset_age)
		elif is_null(cdna1, prot1) or  is_null(cdna2, prot2):
			age_set["one null"].append(onset_age)
		else:
			age_set["none null"].append(onset_age)
	cursor.close()
	db.close()

	fig, axs = plt.subplots(2, 2)
	for n in range(len(age_set_names)):
		name = age_set_names[n]
		axs.flat[n].hist(age_set[name], bins)
		axs.flat[n].legend([f"{name} ({len(age_set[name])})"])

	plt.show()


#########################################
if __name__ == '__main__':
	main()
