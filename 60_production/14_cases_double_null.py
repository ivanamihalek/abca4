#!/usr/bin/python3

import os, re, subprocess
from random import random

from utils.mysql import *
from utils.utils import panic, unpack_progression
import matplotlib.pyplot as plt
from utils.simulation import *


def plot_sim_results_vs_data(all_age, all_va, varid1, varid2, params1, params2, rpe_baseline):

	# simulate
	x, y = sim(params1, params2, rpe_baseline, max_age=50)
	# plot_

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
	qry = "select id, allele_id_1, allele_id_2, onset_age, progression, publication_id, patient_xref_id from cases where  "
	qry += "(progression like '%:%' or (onset_age is not null and onset_age>0)) "
	qry += "and (notes is null or notes not like '%caveat%')"

	all_points_age = []
	all_points_va = []
	patient_ids = {}
	# find the variants corresponding to those cases
	for case_data in hard_landing_search(cursor,qry):
		[case_id, allele_id_1, allele_id_2, onset_age, progression, publication_id, patient_xref_id] = case_data
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
		# print()
		# print(case_id, "onset age:", onset_age)
		# print(progression, publication_id)
		# for ai in [allele_id_1, allele_id_2]:
		# 	print("\t", ai, variants[ai])
		# 	for v in variants[ai]:
		# 		print(f"\t\t variant {v}    {params[v]}")
		if not  publication_id in patient_ids: patient_ids[publication_id] = []
		patient_ids[publication_id].append(patient_xref_id)
		[age, va] = unpack_progression(progression) if progression else [[],[]]
		if onset_age and onset_age>0:
			age = [max(onset_age-1, 0)] + age
			va = [1.0] + va
		all_points_age.append(age)
		all_points_va.append(va)

	for publication_id, patient_xref_ids in patient_ids.items():
		[ref, pubmed] = hard_landing_search(cursor, f"select reference, pubmed from publications where id={publication_id}")[0]
		refpieces = ref.replace("("," ").replace(")"," ").split(" ")
		if len(refpieces)==1: refpieces.append("")
		refstring = "~\\cite{%s%s}:%s" % (refpieces[0].lower(), refpieces[1], ','.join(patient_xref_ids))
		print(refstring, end="; ")
	print()
	plot_sim_results_vs_data(all_points_age, all_points_va, 1, 1, [0.0, 1.0],  [0.0, 1.0], 0.1)

	cursor.close()
	db.close()




#########################################
if __name__ == '__main__':
	main()
