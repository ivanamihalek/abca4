#!/usr/bin/python3

from utils.mysql import *
from utils.abca4_mysql import  *

import matplotlib.pyplot as plt

def variant_info(cursor, allele_id):
	qry = f"select variant_ids from alleles where id={allele_id}"
	variants = []
	for vid in hard_landing_search(cursor, qry)[0][0].strip("-").split("-"):
		qry  = "select p.expression_folding_membrane_incorporation, p.transport_efficiency, v.protein_domain "
		qry += "from variants v, parametrization p "
		qry += f"where p.variant_id=v.id and v.id={vid}"
		ret = hard_landing_search(cursor, qry)
		if len(ret)>2:
			print(f"multiple entries fo variant {vid} (?!)")
			exit()
		variants.append(ret[0])
	return variants


#########################################
def main():

	db, cursor = abca4_connect()
	qry  = "select id, allele_id_1, allele_id_2,  onset_age, progression from cases "
	qry += "where onset_age>0 and (notes is null or notes not like '%caveat%')"
	ages  = {"ECD1":[], "ECD2":[]}
	cases = {"ECD1":[], "ECD2":[]}
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2,  onset_age, progression] = case
		#if onset_age>9: continue
		variants = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			variants.extend(ret)
		if len(variants)!=2: continue
		for ecd in ["ECD1", "ECD2"]:
			if (variants[0][2]==ecd and  variants[1][0]==0.0)  \
						or (variants[1][2]==ecd and  variants[0][0]==0.0):
				ages[ecd].append(onset_age)
				cases[ecd].append(case_id)
	db.close()


	number_of_bins = 75
	bin_size = 1
	bins = []
	for b in range(number_of_bins):
		bins.append(b*bin_size+0.1)

	rows = 2
	cols = 1
	fig, axs = plt.subplots(rows, cols)
	fig.text(0.5, 0.04, 'onset age', ha='center')
	fig.text(0.04, 0.5, 'number of patients', va='center', rotation='vertical')
	n = -1
	for ecd in ["ECD1", "ECD2"]:
		n += 1
		axs.flat[n].hist(ages[ecd], bins)
		axs.flat[n].legend([ecd])
		axs.flat[n].axvline(x=5, color="red")
		axs.flat[n].axvline(x=10, color="red")

	plt.show()

#########################################
if __name__ == '__main__':
	main()


'''
is parametrization consistent with the pehontype produced by variants (are come variants effectively null?)
then run simulation for the initial parametrization and see how far off we are from the reported acuity of vision
does my cire sim represent a single day?
read in our patients to the database

'''