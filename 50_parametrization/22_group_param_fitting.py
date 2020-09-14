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
		variants[vid] = ret[0]
	return variants


def read_in_cases_w_progression():
	db, cursor = abca4_connect()

	# fitting on our own cases on purpose 52 is our set of data
	# take only cases with at least three progression points
	qry  = "select id, allele_id_1, allele_id_2,  progression from cases "
	qry += "where (notes is null or notes not like '%caveat%') "
	qry += "and progression like '%:%:%:%' and publication_id = 52"

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

	cursor.close()
	db.close()

	return [case_variants, prog]

def find_vid_w_exp_parametrization():
	db, cursor = abca4_connect()
	experimental = [str(r[0]) for r in hard_landing_search(cursor, "select variant_id from  parametrization_literature")]
	cursor.close()
	db.close()
	return experimental

#########################################
def find_variant_groups(case_variants):

	# how many cases does a variant belong to?
	variant_cases = {}
	for case_id, variants in case_variants.items():
		for vid in variants.keys():
			if vid not in variant_cases: variant_cases[vid] = []
			variant_cases[vid].append(case_id)

	groups = []
	for case_id, variants in case_variants.items():
		# if both variants belong to only one case (which should hopefully be one and the same)
		# we are not interested
		if sum([len(variant_cases[vid]) for vid in variants.keys()])<3: continue

		# I do not want a case to be joind to a group
		# only through a variant that is nonexpressing - we are not tweaking those
		var_ids = set([vid for vid in variants.keys() if variants[vid][0]>0.001])
		if len(var_ids)==0: continue

		case_belongs_to = []
		for group in groups:
			if len(group.intersection(var_ids))>0:
				case_belongs_to.append(groups.index(group))

		if len(case_belongs_to)==0: # create new_group
			groups.append(var_ids)

		elif len(case_belongs_to)==1: # update existing group
			groups[case_belongs_to[0]].update(var_ids)

		elif len(case_belongs_to)==2: # join two groups
			group1 = groups[case_belongs_to[0]]
			group1.update(var_ids)
			# add the second group to the first
			group2 = groups[case_belongs_to[1]]
			group1.update(group2)
			# and remove the second group
			groups.remove(group2)

		else:
			print(f"case {case_id} belongs to more than 2 groups")
			print(var_ids)
			print(groups)
			exit()

	for group in groups:
		print()
		print("group   ", group)
		print("cases")
		for vid in group:
			print(f"\t {vid}    {variant_cases[vid]}")
	print()
	return groups


def is_in_group(group, case_id, case_variants):
	var_ids = set(list(case_variants[case_id].keys()))
	return len(group.intersection(var_ids))>0



#########################################
from scipy.interpolate import interp1d


def target_function(params, case_ids, case_variants, rpe_baseline, progression):

	d = 0
	norm = 0
	literal_value = 0
	for case_id in case_ids:
		# age,va = last point in progression - the age to which we want to run the simulation
		age, va = [float(number) for number in progression[case_id].split(";")[-1].split(":")]

		# expected va: simulate to  nearest rounded age
		vid = list(case_variants[case_id].keys())
		x, y = sim(params[vid[0]], params[vid[1]], rpe_baseline[case_id], int(round(age))+5)
		interpol = interp1d(x, y["rpe"])
		for agepoint in progression[case_id].split(";"):
			age, va = [float(p) for p in agepoint.split(":")]
			expected_va = interpol([age])[0]
			d += (va - expected_va)**2
			norm += 1
			literal_value += expected_va - va

	d = sqrt(d/norm) if norm>0 else 0
	literal_value = literal_value/norm if norm>0 else 0
	parametrization_descriptor = {"rmsd": d, "sign": literal_value}
	return parametrization_descriptor


########################################
# the probabilty that we should tweak one particular characteristic value of a variant
def cutoff_prob_by_annotation(notes):
	cutoff_prob = 0.5
	if 'expression' in notes: # it is more likely that we should tweak the expression
		if 'no expression' in notes:
			cutoff_prob = 0.9
		elif 'severe expression' in notes:
			cutoff_prob = 0.8
		elif 'strong expression' in notes:  # means strongly affected, not robust expression
			cutoff_prob = 0.7
		elif 'mild expression' in notes:
			cutoff_prob = 0.6

	elif 'transport' in notes:
		if 'no transport' in notes:
			cutoff_prob = 0.1
		elif 'severe transport' in notes:
			cutoff_prob = 0.2
		elif 'strong transport' in notes: # means strongly affected, not robust expression
			cutoff_prob = 0.3
		elif 'mild transport' in notes:
			cutoff_prob = 0.4

	return cutoff_prob

def tweak_variant_params(params, character, vid, is_experimental):
	step_size = 0.005 if is_experimental else 0.01
	[e_prev, t_prev] = params[vid]
	[cdna, prot, notes] = character[vid]
	[e,t] = [e_prev, t_prev]
	# if dist_per_var<0:
	die = random()
	if die < 0.5: #  tune up one of the values
		# splice mutation does not change transport properties
		if 'splice' in prot:
			e = e_prev + step_size
		else:
			cutoff_prob = cutoff_prob_by_annotation(notes)
			if random()<cutoff_prob:
				e = e_prev + step_size
			else:
				t = t_prev + step_size
		if e>1.0: e=1.0
		if t>1.0: t=1.0

	else:  # tune down one of the values
		if 'splice' in prot:
			e = e_prev - step_size
		else:
			cutoff_prob = cutoff_prob_by_annotation(notes)
			if random()<cutoff_prob:
				e = e_prev - step_size
			else:
				t = t_prev - step_size
		if e<0.0: e=0.0
		if t<0.0: t=0.0

	params[vid] = [e,t]

	return


########################################
def tweak_baseline(rpe_baseline, case_id):
	# pick case at random and tweak its baseline
	if 0 < rpe_baseline[case_id] <= 0.2:
		if random() < 0.5:
			rpe_baseline[case_id] = max(0.0, rpe_baseline[case_id]-0.001)
		else:
			rpe_baseline[case_id] = min(0.2, rpe_baseline[case_id]+0.001)

	elif rpe_baseline[case_id]>0.2: #try to decrease
		rpe_baseline[case_id] = max(0.0, rpe_baseline[case_id]-0.001)

	else:
		rpe_baseline[case_id] = min(0.2, rpe_baseline[case_id]+0.001)


##########
def accept_or_reject_parameter_change(parametrization, case_ids, case_variants, progression, identifiers,  prev):
	[params, rpe_baseline, T] = parametrization
	[vid, case_id] = identifiers
	[min_dist, min_params, min_rpe_baseline, e_prev, t_prev, rpe_baseline_prev, dist_prev] = prev
	param_descr = target_function(params, case_ids, case_variants, rpe_baseline, progression)
	dist = param_descr["rmsd"]
	# print(" %7.2f        %7.2f   %7.2f  %7.2f         %7.2f   %7.2f   %7.2f         %7.3f   %7.3f   " %
	#       (param_descr_prev["sign"],  e_prev, t_prev, rpe_baseline_prev, e, t, rpe_baseline, dist_prev, dist))
	if min_dist>dist:
		min_dist = dist
		min_params = copy.deepcopy(params)
		min_rpe_baseline = rpe_baseline

	if dist<=dist_prev:
		# accept
		dist_prev = dist
	# elif random()< exp(-(dist-dist_prev)/T):
	# 	# print("accepted %.2f --> %.2f      %.3f" % (dist_prev, dist, exp(-(dist-dist_prev)/T)))
	# 	# accept
	# 	dist_prev = dist
	else:
		# backpedal
		if vid: params[vid] = [e_prev, t_prev]
		if case_id: rpe_baseline[case_id] = rpe_baseline_prev

	return [min_dist, min_params, min_rpe_baseline, dist_prev]


##########
def optimization_loop(orig_params, case_ids, case_variants, orig_rpe_baseline, progression, character,  experimental, verbose=True):

	# variants which we know are null at the start stay null
	# some variants may become null during simulation, and
	# we want to distiguish the two
	# (but we do want the sitance sum using these fixed variants)
	fixed = {}
	for vid, prms in orig_params.items():
		fixed[vid] = prms[0]<0.0001

	params = copy.deepcopy(orig_params)
	rpe_baseline = copy.deepcopy(orig_rpe_baseline)

	T = 0.01
	min_dist = 10
	min_params = copy.deepcopy(orig_params)
	min_rpe_baseline = copy.deepcopy(orig_rpe_baseline)

	for pass_number in range(1,500):

		param_descr_prev = target_function(params, case_ids, case_variants, rpe_baseline, progression)
		dist_prev = param_descr_prev["rmsd"]

		for vid, prms in params.items():
			if fixed[vid]: continue # we do not change the null variants
			#if vid in experimental: continue
			[e_prev, t_prev] = params[vid]
			tweak_variant_params(params, character, vid, vid in experimental)
			prev = [min_dist, min_params, min_rpe_baseline, e_prev, t_prev, None, dist_prev]
			ret = accept_or_reject_parameter_change([params, rpe_baseline, T], case_ids, case_variants, progression, [vid, None], prev)
			[min_dist, min_params, min_rpe_baseline, dist_prev] = ret

		for case_id in case_ids:
			rpe_baseline_prev = rpe_baseline[case_id]
			tweak_baseline(rpe_baseline, case_id)
			prev = [min_dist, min_params, min_rpe_baseline, None, None, rpe_baseline_prev, dist_prev]
			ret = accept_or_reject_parameter_change([params, rpe_baseline, T], case_ids, case_variants, progression, [None, case_id], prev)
			[min_dist, min_params, min_rpe_baseline, dist_prev] = ret

	if verbose: print("min dist found:  %7.3f  " % min_dist)
	return min_params, min_rpe_baseline

#########################################
def fit_to_variant_group(group, case_ids, case_variants, progression, experimental):
	rpe_baseline_init = 0.1
	print(group)

	# the params shoud be per variant
	params = {}
	character = {}
	rpe_baseline = {}
	for case_id in case_ids:
		rpe_baseline[case_id] = rpe_baseline_init
		for var_id, var_descr in case_variants[case_id].items():
			if var_id in params: continue
			params[var_id] = var_descr[:2]
			character[var_id] = var_descr[2:5]

	return optimization_loop(params, case_ids, case_variants, rpe_baseline, progression, character, experimental)

#########################################
def store_variant_params(new_variant_params, group_id):
	db, cursor = abca4_connect()
	for vid, newprms in new_variant_params.items():
		[enew, tnew] = newprms
		store_or_update(cursor, 'parametrization_adjusted', fixed_fields={"variant_id":vid},
		                update_fields={"expression_folding_membrane_incorporation": enew,
		                               "transport_efficiency": tnew, "notes": f"group {group_id}"})
	cursor.close()
	db.close()

#########################################
def store_baseline(new_baseline):
	db, cursor = abca4_connect()
	for case_id, baseline in new_baseline.items():
		store_or_update(cursor, 'case_parametrization', fixed_fields={"case_id":case_id},
		                update_fields={"baseline":baseline})
	cursor.close()
	db.close()


#########################################
def main():

	[case_variants, progression] = read_in_cases_w_progression()
	groups = find_variant_groups(case_variants)
	experimental = find_vid_w_exp_parametrization()

	case_ids_belonging_to_group = {}
	for g in range(len(groups)):
		group = groups[g]
		if g not in case_ids_belonging_to_group: case_ids_belonging_to_group[g] = []
		for case_id, vars in case_variants.items():
			if not is_in_group(group, case_id, case_variants): continue
			case_ids_belonging_to_group[g].append(case_id)

	for g, case_ids in case_ids_belonging_to_group.items():
		group = groups[g]
		[new_variant_params, new_rpe_baseline]  = fit_to_variant_group(group, case_ids, case_variants, progression, experimental)
		store_variant_params(new_variant_params, g+1)
		store_baseline(new_rpe_baseline)


#########################################
if __name__ == '__main__':
	main()

