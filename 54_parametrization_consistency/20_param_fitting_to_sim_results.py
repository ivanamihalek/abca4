#!/usr/bin/python3
import copy
from math import sqrt
from random import random

from utils.abca4_mysql import  *

from utils.simulation import progression_sim
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


def read_in_values():
	db, cursor = abca4_connect()
	qry  = "select id, allele_id_1, allele_id_2,  progression from cases "
	qry += "where (notes is null or notes not like '%caveat%') and progression is not null and progression!=''"
	# we don't want to parametrize on our own cases, but they are already filtered out because they do not contain the onset age
	# (the onset age is set to -1)

	# fitting on our own cases on purpose
	# qry  = "select id, allele_id_1, allele_id_2,  progression from cases "
	# qry += "where (notes is null or notes not like '%caveat%') and publication_id = 52"

	prog = {}
	case_variants = {}
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2, progression ] = case
		#if len(progression.split(";"))<3: continue
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

def sim(var1_param, var2_param, max_age):
	# the values for one case/patient
	# alpha is the abilty to fold and incorporate
	# "fraction" refers to the fraction of the wild-type capability
	expression_fraction = [var1_param[0], var2_param[0]]
	transport_fraction = [var1_param[1], var2_param[1]]
	x, y  = progression_sim(max_age, expression_fraction, transport_fraction, verbose=False)
	return x, y


#########################################
def target_function_0(params, case_variants, progression):
	d = 0
	norm = 0
	for case_id, vars in case_variants.items():
		# age,va = fist point in progression
		age, va = [float(number) for number in progression[case_id].split(";")[0].split(":")]
		# expected va: simulate to  nearest rounded age
		vid = list(vars.keys())
		if len(vid)!=2: panic("at the moemnt I can handle one and only one var per allele")
		x, y = sim(params[vid[0]], params[vid[1]], int(round(age))+1)
		last_age_in_sim = x[-1]
		expected_va = y["rpe"][-1]
		# print("\n", case_id, age, va)
		# print(f"expected va at age {last_age_in_sim}: %.2f"%expected_va)
		d += (va - expected_va)**2
		norm += 1
	return sqrt(d/norm) if norm>0 else 0

#########################################
def target_function(params, case_variants, progression, character=None):
	d = 0
	norm = 0
	# collect the top offenders as we go
	rmsd_distances_per_variant = {}
	distances_per_variant = {}
	for case_id, vars in case_variants.items():
		# age,va = fist point in progression
		age, va = [float(number) for number in progression[case_id].split(";")[0].split(":")]
		# expected va: simulate to  nearest rounded age
		vid = list(vars.keys())
		if len(vid)!=2: panic("at the moemnt I can handle one and only one var per allele")
		x, y = sim(params[vid[0]], params[vid[1]], int(round(age))+1)
		last_age_in_sim = x[-1]
		expected_va = y["rpe"][-1]
		# print("\n", case_id, age, va)
		# print(f"expected va at age {last_age_in_sim}: %.2f"%expected_va)
		d += (va - expected_va)**2
		norm += 1

		for v in vid:
			if not v in distances_per_variant:
				distances_per_variant[v] = 0
				rmsd_distances_per_variant[v] = 0
			rmsd_distances_per_variant[v] += (expected_va - va)**2
			distances_per_variant[v] += (expected_va - va)

	d = sqrt(d/norm) if norm>0 else 0
	variants_sorted_by_distance = sorted(rmsd_distances_per_variant.keys(),
	                                     key = lambda v: rmsd_distances_per_variant[v], reverse=True)

	if character:
		for v in variants_sorted_by_distance[:20]:
			print(" %3d    %7.3f    %7.3f   " %(int(v), rmsd_distances_per_variant[v], distances_per_variant[v]),
			      "   ", params[v],  "   ", character[v])

	parametrization_descriptor = {"rmsd": d, "vars_by_dist": variants_sorted_by_distance,
	                              "rmsd_per_var": rmsd_distances_per_variant,
	                              "dist_per_var": distances_per_variant}
	return parametrization_descriptor


########################################
def tweak(e_prev, t_prev, dist_per_var, character):
	[cdna, prot, notes] = character
	[e,t] = [e_prev, t_prev]
	if dist_per_var<0:
		# we need to increase the (presumed) performance
		if 'splice' in prot:
			e = e_prev + 0.2
		else:
			if random()<0.5:
				e = e_prev + 0.1
			else:
				t = t_prev+ 0.1

		if e>1.0: e=1.0
		if t>1.0: t=1.0


	else:
		if 'splice' in prot:
			e = e_prev - 0.2
		else:
			if random()<0.5:
				e = e_prev - 0.1
			else:
				t = t_prev - 0.1

		if e<0.0: e=0.0
		if t<0.0: t=0.0


	return [e,t]

##########
def optimization_loop(orig_params, case_variants, progression, character):
	case_id = list(case_variants.keys())
	number_of_cases = len(case_id)

	# variants which we know are null at the start stay null
	# some variants may become null during simulation, and
	# we want to distiguish the two
	# (but we do want the sitance sum using these fixed variants)
	fixed = {}
	for vid, prms in orig_params.items():
		fixed[vid] = prms[0]<0.0001

	params = copy.deepcopy(orig_params)

	T = 1
	for pass_number in range(1,30):
		param_descr_prev = target_function(params, case_variants, progression)
		dist_prev = param_descr_prev["rmsd"]
		print("%2d     %.2f"%(pass_number, dist_prev))

		ct = 0
		top_offenders = min(pass_number*10, 100)
		#for vid, prms in params.items():
		for vid in param_descr_prev["vars_by_dist"][:top_offenders]:
			if fixed[vid]: continue # we do not change the null variants
			ct += 1
			if not ct%50: print(f"{ct} out of {len(params)}")
			[e_prev, t_prev] = params[vid]
			[e, t] = tweak(e_prev, t_prev, param_descr_prev["dist_per_var"][vid], character[vid])
			params[vid] = [e,t]
			param_descr = target_function(params, case_variants, progression)
			dist = param_descr["rmsd"]
			print(" %7.2f   %7.2f   %7.2f    %7.2f   %7.2f      %7.3f   %7.3f   " %
			      (param_descr_prev["dist_per_var"][vid],  e_prev, t_prev, e, t, dist_prev, dist), character[vid])
			if dist<dist_prev:
				# accept
				dist_prev = dist
			else:
				# backpedal

				params[vid] = [e_prev,t_prev]

		print("%2d     %.2f\n"%(pass_number, dist_prev))

	return params


def store(parameters, new_variant_params, character):
	db, cursor = abca4_connect()
	for vid, newprms in new_variant_params.items():
		[eold, told] = parameters[vid]
		[enew, tnew] = newprms
		print( "%3d  %7.2f   %7.2f   %7.2f   %7.2f  " % (int(vid), eold, told, enew, tnew), character[vid])
		store_or_update(cursor, 'parametrization_adjusted', fixed_fields={"variant_id":vid},
		                update_fields={"expression_folding_membrane_incorporation":enew,
		                               "transport_efficiency":tnew})
	cursor.close()
	db.close()

#########################################
def main():

	[case_variants, progression] = read_in_values()
	parameters = {}
	character = {}
	for case_id, vars in case_variants.items():
		for vid, variant_info in vars.items():
			parameters[vid]  = variant_info[:2]
			character[vid] = variant_info[2:5]

	print(f"number of vars {len(parameters)}")

	param_descr  = target_function(parameters, case_variants, progression, character)
	print("intial d = %.2f"% param_descr["rmsd"])

	new_variant_params = optimization_loop(parameters, case_variants, progression, character)

	store(parameters, new_variant_params, character)

#########################################
if __name__ == '__main__':
	main()

