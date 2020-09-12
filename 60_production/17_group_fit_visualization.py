#!/usr/bin/python3
import copy
from math import sqrt, exp
from random import random, choice

from utils.abca4_mysql import  *
from utils.utils import unpack_progression
from utils.simulation import *
import matplotlib.pyplot as plt


def panic(msg):
	print(msg)
	exit(1)


def variant_info(cursor, allele_id):
	qry = f"select variant_ids from alleles where id={allele_id}"
	variants = {}
	for vid in hard_landing_search(cursor, qry)[0][0].strip("-").split("-"):
		qry  = "select p.expression_folding_membrane_incorporation, p.transport_efficiency,   "
		qry += "v.cdna, v.protein, p.notes "
		qry += "from variants v, parametrization p "
		qry += f"where p.variant_id=v.id and v.id={vid}"
		ret = hard_landing_search(cursor, qry)
		if len(ret)>2:
			print(f"multiple entries fo variant {vid} (?!)")
			exit()
		variants[int(vid)] = ret[0]
	return variants

####################################
def read_progression(cursor):
	# fitting on our own cases on purpose 52 is our set of data
	qry  = "select id, allele_id_1, allele_id_2,  progression from cases "
	qry += "where (notes is null or notes not like '%caveat%') and publication_id = 52"

	prog = {}
	case_variants = {}
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2, progression ] = case
		if len(progression.split(";"))<3: continue
		variants = {}
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			variants.update(ret)
		if len(variants)!=2: continue
		if not case_id in case_variants: case_variants[case_id] = {}
		case_variants[case_id].update(variants)
		prog[case_id] = progression

	return [case_variants, prog]

####################################
def read_baseline(cursor):
	baseline = {}
	old_baseline = {}
	for case_id, bsln in hard_landing_search(cursor, "select case_id, baseline from case_parametrization"):
		baseline[case_id] = bsln
		old_baseline[case_id] = 0.1

	return [baseline, old_baseline]

###########
def read_groups(cursor):
	group_variants = {}
	params = {}
	old_params = {}
	qry = "select * from parametrization_adjusted where notes like 'group%'"
	for [prm_id, variant_id, expression, transport, notes] in  hard_landing_search(cursor, qry):
		group = int(notes.split()[-1])
		if not group in group_variants: group_variants[group] = []
		group_variants[group].append(variant_id)
		params[variant_id] = [expression, transport]
		qry  = "select expression_folding_membrane_incorporation, transport_efficiency "
		qry += f"from parametrization where variant_id={variant_id}"
		old_params[variant_id] = hard_landing_search(cursor,qry)[0]

	return [group_variants, params, old_params]

#########################################
def report(parameters, rpe_baseline, new_variant_params, new_rpe_baseline, varid, case_id):
	print("\n========================")
	print(f"case id {case_id}")

	print(f"{varid[0]}:  %.2f  %.2f\n{varid[1]}:  %.2f  %.2f\nrpe_baseline:   %.2f" % \
	      (parameters[varid[0]][0], parameters[varid[0]][1],
	       parameters[varid[1]][0], parameters[varid[1]][1], rpe_baseline[case_id]))
	print("--------------------------")

	print(f"{varid[0]}:  %.2f  %.2f\n{varid[1]}:  %.2f  %.2f\nrpe_baseline:   %.2f" % \
	      (new_variant_params[varid[0]][0], new_variant_params[varid[0]][1],
	       new_variant_params[varid[1]][0], new_variant_params[varid[1]][1], new_rpe_baseline[case_id]))
	print()


def multiplot_sim_results(axs, index, progression, params_old, params_new,
                          rpe_baseline, rpe_baseline_new, varid, case_id):

	age_va, va = unpack_progression(progression[case_id])

	# simulate old
	params1 = params_old[varid[0]]
	params2 = params_old[varid[1]]
	x_old, y_old = sim(params1, params2, rpe_baseline[case_id], max_age=50)

	# simulate new
	params1 = params_new[varid[0]]
	params2 = params_new[varid[1]]
	x, y = sim(params1, params2, rpe_baseline_new[case_id], max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	title = f"{varid[0]}:  %.2f  %.2f\n{varid[1]}:  %.2f  %.2f" % (params1[0], params1[1], params2[0], params2[1])

	axs.flat[index].set_ylim(-0.1,1.1)
	axs.flat[index].set_title(title)
	axs.flat[index].plot(x, y["rpe"])
	axs.flat[index].plot(x_old, y_old["rpe"], linestyle="dashed")
	axs.flat[index].scatter(age_va, va, color='red')

	return

def visualization(case_ids, case_variants, progression, params, rpe_baseline_new, params_old, rpe_baseline_old):

	if len(case_ids)<3:
		rows = 1
		cols = 2
	elif len(case_ids)<5:
		rows = 2
		cols = 2
	else:
		rows = 2
		cols = 4

	fig, axs = plt.subplots(rows, cols)
	index = -1
	for case_id in case_ids:
		varids = list(case_variants[case_id].keys())
		report(params_old, rpe_baseline_old, params, rpe_baseline_new, varids, case_id)
		index += 1
		multiplot_sim_results(axs, index, progression, params_old, params, rpe_baseline_old, rpe_baseline_new, varids, case_id)

	plt.show()



#########################################
def main():

	db, cursor = abca4_connect()

	[case_variants, progression] = read_progression(cursor)
	[groups, params, old_params] = read_groups(cursor)
	[rpe_baseline, old_rpe_baseline] = read_baseline(cursor)
	cursor.close()
	db.close()

	case_ids_belonging_to_group = {}
	for g, group_varids in groups.items():
		if g not in case_ids_belonging_to_group: case_ids_belonging_to_group[g] = []
		for case_id, vars in case_variants.items():
			for varid in vars:
				if varid in group_varids:
					case_ids_belonging_to_group[g].append(case_id)
					break

	for g, varids in groups.items():
		case_ids = case_ids_belonging_to_group[g]
		if len(case_ids )<2: continue
		visualization(case_ids, case_variants, progression, params, rpe_baseline, old_params, old_rpe_baseline)

#########################################
if __name__ == '__main__':
	main()

