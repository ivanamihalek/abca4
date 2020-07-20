#!/usr/bin/python3 -u


from utils.utils import  *
from utils.abca4_mysql import  *

import matplotlib.pyplot as plt

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
		[cdna1, prot1, freqs1, homozygs1, region1, cons1] = variants_from_allele(cursor, allele_id_1)
		[cdna2, prot2, freqs2, homozygs2, region2, cons2] = variants_from_allele(cursor, allele_id_2)
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



#########################################
if __name__ == '__main__':
	main()
