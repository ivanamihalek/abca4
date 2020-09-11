#!/usr/bin/python3

import os, re, subprocess
from random import random

from utils.mysql import *
from utils.utils import panic, unpack_progression
import matplotlib.pyplot as plt
from utils.simulation import *


named_colors = ['palegreen', 'lime', 'green', 'salmon', 'indianred', 'darkred', 'sienna',
                'tomato', 'rosybrown', 'darkorange', 'orange', 'khaki', 'gold', 'yellow', 'greenyellow',
                'red', 'lightcoral', 'firebrick']

def plot_sim_results_vs_data(all_age, all_va, varid1, varid2, params1, params2, rpe_baseline):

	# simulate
	x, y = sim(params1, params2, rpe_baseline, max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	title = f"{varid1}:  %.2f  %.2f\n{varid2}:  %.2f  %.2f" % (params1[0], params1[1], params2[0], params2[1])

	rows = 4
	cols = 7

	fig, axs = plt.subplots(rows, cols)

	for i in range(len(all_age)):
		print(i, all_age[i], all_va[i])
		axs.flat[i].set_ylim(-0.1,1.1)
		axs.flat[i].set_xticklabels([])
		axs.flat[i].set_yticklabels([])

		axs.flat[i].plot(x, y["rpe"])
		axs.flat[i].scatter(all_age[i],  all_va[i], color='red')

	plt.show()
	return

#########################################
def main():

	ter_only = False
	db, cursor = abca4_connect()

	# find cases with at least three progression points
	qry = "select id, allele_id_1, allele_id_2, onset_age, progression, publication_id from cases where  "
	qry += "(progression like '%:%' or (onset_age is not null and onset_age>0)) "
	qry += "and (notes is null or notes not like '%caveat%')"

	all_points_age = []
	all_points_va = []
	# find the variants corresponding to those cases
	for [case_id, allele_id_1, allele_id_2, onset_age, progression, publication_id] in hard_landing_search(cursor,qry):
		params = {}
		variants = {}
		for ai in [allele_id_1, allele_id_2]:
			variants[ai] = hard_landing_search(cursor, f"select variant_ids from alleles where id={ai}")[0][0].strip("-").split("-")
			for v in variants[ai]:
				# is this null by any chance?
				if ter_only:
					qry = f"select * from variants where id={v} and protein like '%Ter%'"
				else:
					qry  = "select * from variants v, parametrization p "
					#qry += f"where v.id={v} and v.protein not like '%Ter%' "
					qry += f"where v.id={v} and (v.protein  like '%Ter%' or  v.protein  like 'splice') "
					qry += f"and p.variant_id={v} and p.expression_folding_membrane_incorporation<0.001"
				ret = error_intolerant_search(cursor, qry)
				if ret:
					# if len(ret)>1: panic([f"multiple parametrization for varid {v}"])
					# [prm_id, var_id, e, t, notes] = ret[0]
					# if e>0.001: continue # this is expressing
					# params[v] = [prm_id, var_id, e, t, 0]
					params[v] = [1, 1, 0.0, 1.0, 0]

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

		[age, va] = unpack_progression(progression) if progression else [[],[]]
		if onset_age and onset_age>0:
			age = [max(onset_age-1, 0)] + age
			va = [1.0] + va
		all_points_age.append(age)
		all_points_va.append(va)

	plot_sim_results_vs_data(all_points_age, all_points_va, 1, 1, [0.0, 1.0],  [0.0, 1.0], 0.1)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
