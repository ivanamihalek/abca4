#!/usr/bin/python3

import os, re, subprocess

from utils.mysql import *
from utils.utils import panic, unpack_progression
import matplotlib.pyplot as plt
from utils.simulation import *


def plot_sim_results_vs_data(age, va, varid1, varid2, params1, params2, rpe_baseline):

	# simulate
	x, y = sim(params1, params2, rpe_baseline, max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	title = f"{varid1}:  %.2f  %.2f\n{varid2}:  %.2f  %.2f" % (params1[0], params1[1], params2[0], params2[1])

	fig, axs = plt.subplots()
	axs.set_ylim(-0.1,1.1)
	axs.set_title(title)
	axs.plot(x, y["rpe"])
	axs.scatter(age, va, color='red')

	plt.show()
	return

#########################################
def main():


	db, cursor = abca4_connect()

	# find cases with at least three progression points
	qry = "select id, allele_id_1, allele_id_2, onset_age, progression, publication_id from cases where  "
	qry += "(progression like '%:%:%:%' or (progression like '%:%:%' and onset_age is not null and onset_age>0)) "
	qry += "and (notes is null or notes not like '%caveat%')"

	# find the variants corresponding to those cases
	for [case_id, allele_id_1, allele_id_2, onset_age, progression, publication_id] in hard_landing_search(cursor,qry):
		params = {}
		variants = {}
		for ai in [allele_id_1, allele_id_2]:
			variants[ai] = hard_landing_search(cursor, f"select variant_ids from alleles where id={ai}")[0][0].strip("-").split("-")
			for v in variants[ai]:
				ret = error_intolerant_search(cursor, f"select * from parametrization_adjusted where variant_id={v}")
				if False and ret:
					if len(ret)>1: panic([f"multiple parametrizations for varid {v}"])
					params[v] = ret[0]
					continue
				else:
					ret = error_intolerant_search(cursor, f"select * from parametrization_literature where variant_id={v}")
					if not ret: continue
					if len(ret)>1: panic([f"multiple parametrizations for varid {v}"])
					params[v] = ret[0]

		if len(params)!=(len(variants[allele_id_1])+len(variants[allele_id_2])): continue

		# keep only if all variants have experimental  support
		# i am still not ready to deal with multiple variants per allele
		if len(variants[allele_id_1])>1 or len(variants[allele_id_2])>1: continue
		print()
		print(case_id, "onset age:", onset_age)
		print(progression, publication_id)
		for ai in [allele_id_1, allele_id_2]:
			print("\t", ai, variants[ai])
			for v in variants[ai]:
				print(f"\t\t variant {v}    {params[v]}")
		varid1 = variants[allele_id_1][0]
		varid2 = variants[allele_id_2][0]
		params1 = params[varid1][2:4]
		params2 = params[varid2][2:4]
		print(varid1, params1)
		print(varid2, params2)
		age, va = unpack_progression(progression)
		if onset_age and onset_age>0:
			age = [max(onset_age-1, 0)] + age
			va = [1.0] + va
		# plot_sim_results_vs_data(age, va, varid1, varid2, params1, params2, rpe_baseline)
		plot_sim_results_vs_data(age, va, varid1, varid2, params1, params2, 0.1)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
