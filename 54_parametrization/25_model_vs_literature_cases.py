#!/usr/bin/python3
import copy
from math import sqrt, exp
from random import random, choice

from utils.abca4_mysql import *
from utils.utils import unpack_progression

from utils.simulation import *
import matplotlib.pyplot as plt

#########################################
def var_allele_mapping(cursor, params):
	var2allele = {}
	allele2var = {}
	for vid in params.keys():
		qry = f"select id from alleles where variant_ids='-{vid}-'"
		var2allele[vid] = hard_landing_search(cursor, qry)[0][0]
		allele2var[var2allele[vid]] = vid
	return var2allele, allele2var


#########################################
def var_characterization(cursor, params):
	characterization = {}
	for vid in params.keys():
		qry = f"select notes from parametrization where variant_id={vid}"
		characterization[vid] = hard_landing_search(cursor, qry)[0][0]

	return characterization


def literature_cases(cursor, variant_ids,this_publication_id):
	alleles = {}
	progression = {}
	allele_id_string = ",".join([str(aid) for aid in variant_ids])
	qry  =  "select id, allele_id_1, allele_id_2, onset_age, progression from cases "
	qry += f"where allele_id_1 in ({allele_id_string}) and allele_id_2 in ({allele_id_string}) "
	qry += f"and publication_id!={this_publication_id}"
	for line in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2, onset_age, progr_string] = line
		age, va = unpack_progression(progr_string)
		if onset_age and onset_age>0:
			# age = [max(onset_age-1,0)] + age
			# va  = [1.0] + va
			age = [onset_age] + age
			va  = [0.9] + va
		alleles[case_id] = [allele_id_1, allele_id_2]
		progression[case_id] = [age, va]

	return alleles, progression


#########################################
def plot_sim_results_vs_data(age, va, params1, params2, rpe_baseline):

	# simulate
	x, y = sim(params1, params2, rpe_baseline, max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	title = f"a1:  %.2f  %.2f\na2:  %.2f  %.2f" % (params1[0], params1[1], params2[0], params2[1])

	rows = 1
	cols = 1
	fig, axs = plt.subplots(rows, cols)
	axs.set_ylim(-0.1,1.1) # axs[0].set ... if number of subplots>1
	axs.set_title(title)
	axs.plot(x, y["rpe"])
	axs.scatter(age, va, color='red')

	plt.show()
	return

#########################################
def main():
	this_publication_id = 52
	group = 2
	db, cursor = abca4_connect()

	# variants in the group
	qry = "select  variant_id, expression_folding_membrane_incorporation, transport_efficiency "
	qry += f"from parametrization_adjusted where notes like 'group {group}'"
	params = dict([(ret[0], [ret[1], ret[2]]) for ret in hard_landing_search(cursor, qry)])

	# variant to allele mapping
	var2allele, allele2var = var_allele_mapping(cursor, params)
	characterization = var_characterization(cursor, params)

	# literature cases that have these alleles
	alleles, progression = literature_cases(cursor, list(var2allele.keys()), this_publication_id)
	for case_id, [a1, a2] in alleles.items():
		[v1,v2] = [allele2var[a] for a in [a1, a2]]
		print(case_id, a1, a2, v1, v2, progression[case_id])
		print(v1, characterization[v1])
		print(v2, characterization[v2])
		[age, va] =  progression[case_id]
		plot_sim_results_vs_data(age, va, params[v1], params[v2], 0.1)

#########################################
if __name__ == '__main__':
	main()

