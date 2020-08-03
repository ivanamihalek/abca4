#!/usr/bin/python3

from utils.mysql import *
from utils.abca4_mysql import  *
from math import sqrt
from random import random
import copy

def variant_info(cursor, allele_id):
	qry = f"select variant_ids from alleles where id={allele_id}"
	variants = []
	for vid in hard_landing_search(cursor, qry)[0][0].strip("-").split("-"):
		qry  = "select v.id, p.expression_folding_membrane_incorporation, p.transport_efficiency "
		qry += "from variants v, parametrization p "
		qry += f"where p.variant_id=v.id and v.id={vid}"
		ret = hard_landing_search(cursor, qry)
		if len(ret)>2:
			print(f"multiple entries fo variant {vid} (?!)")
			exit()
		variants.append(ret[0])
	return variants

def round_to_4(floatnum):
	return int(round(floatnum*4, 0))


def read_in_values():
	db, cursor = abca4_connect()
	qry  = "select id, allele_id_1, allele_id_2,  onset_age from cases "
	qry += "where onset_age>0 and (notes is null or notes not like '%caveat%')"
	age = {}
	case_variants = {}
	variant_parameters = {}
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2,  onset_age] = case
		if onset_age<0: continue
		variants = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			variants.extend(ret)
		if len(variants)!=2: continue

		case_variants[case_id] = []
		for [vid, expr, transp] in variants:
			# variant_parameters[vid] = [round_to_4(expr), round_to_4(transp)]
			variant_parameters[vid] = f"{round_to_4(expr)}_{round_to_4(transp)}"
			case_variants[case_id].append(vid)

		age[case_id] = onset_age


	db.close()

	# for case_id, varids in case_variants.items():
	# 	print(case_id, age[case_id], varids)
	# print()
	# for varid, params in variant_parameters.items():
	# 	print(varid, params)

	return [case_variants, age, variant_parameters]

#########################################
def same_params(variants1, variants2, variant_parameters):
	params1 = set([variant_parameters[v] for v in variants1])
	params2 = set([variant_parameters[v] for v in variants2])
	return  params1 == params2

########
def rmsd_age_distance(case_variants, age, variant_parameters, case_id, number_of_cases):
	dist = 0
	norm = 0
	for i in range(number_of_cases):
		variants1 = case_variants[case_id[i]]
		a1 = age[case_id[i]]
		for j in range(i+1, number_of_cases):
			variants2 = case_variants[case_id[j]]
			if same_params(variants1, variants2, variant_parameters):
				a2 = age[case_id[j]]
				dist += (a1-a2)**2
				norm += 1
	dist = sqrt(dist/norm)
	return dist

def rmsd_distance_any(inlist):
	dist = 0
	norm = 0
	n = len(inlist)
	for i in range(n):
		a1 = inlist[i]
		for j in range(i+1, n):
			a2 = inlist[j]
			dist += (a1-a2)**2
			norm += 1
	dist = sqrt(dist/norm) if norm else 0
	return dist


def tweak(zero2four):
	if zero2four==0:
		return 1
	if zero2four==4:
		return 3
	if random()<0.5:
		return zero2four+1
	return zero2four-1


##########
def optimization_loop(case_variants, age, orig_variant_parameters):
	case_id = list(case_variants.keys())
	number_of_cases = len(case_id)

	# variants which we know are null at the start stay null
	# some variants may become null during simulation, and
	# we want to distiguish the two
	# (but we do want the sitance sum using these fixed variants)
	fixed = {}
	for v, prms in orig_variant_parameters.items():
		fixed[v] = "0_" in prms

	variant_parameters = copy.deepcopy(orig_variant_parameters)

	# pull the original nulls to the new set
	dist_prev = rmsd_age_distance(case_variants, age, variant_parameters, case_id, number_of_cases)
	print("%2.0f"%dist_prev)
	T = 1
	for pass_number in range(2):
		for v, prms in variant_parameters.items():
			if fixed[v]: continue # we do not change the null variants
			[e_prev, t_prev] =[int(i) for i in prms.split("_")]
			if random()<0.5:
				e = tweak(e_prev)
				t = t_prev
			else:
				e = e_prev
				t = tweak(t_prev)
			variant_parameters[v] = f"{e}_{t}"
			dist = rmsd_age_distance(case_variants, age, variant_parameters, case_id, number_of_cases)
			if dist<dist_prev:
				# accept
				dist_prev = dist
			else:
				# backpedal
				variant_parameters[v] = f"{e_prev}_{t_prev}"

		print("%2d     %.2f"%(pass_number+1, dist_prev))

	return variant_parameters

##############
def plot_distros(case_variants, age, variant_parameters):
	age_list = {}
	for case_id, vars in case_variants.items():

		label = " | ".join(sorted([variant_parameters[v] for v in vars]))
		if not label in age_list: age_list[label] = []
		age_list[label].append(age[case_id])

	labels_sorted = sorted(age_list.keys())
	for label in labels_sorted:
		agelist  = age_list[label]
		dist = rmsd_distance_any(agelist)
		print(label, "       %2.0f     %2.0f"%(sum(agelist)/len(agelist), dist),  "      ", sorted(agelist))
	print(f"number of labels {len(age_list)}")

#########################################
def main():

	[case_variants, age, variant_parameters] = read_in_values()
	plot_distros(case_variants, age, variant_parameters)
	exit()

	new_variant_params = optimization_loop(case_variants, age, variant_parameters)

	plot_distros(case_variants, age, new_variant_params)

	# for v, prms in variant_parameters.items():
	# 	if prms != new_variant_params[v]:
	# 		print(v, prms, new_variant_params[v])



#########################################
if __name__ == '__main__':
	main()

