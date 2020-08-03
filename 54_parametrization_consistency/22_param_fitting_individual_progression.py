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

	# fitting on our own cases on purpose
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
from scipy.interpolate import interp1d

def target_function(params, variants, progression, character=None):

	# collect the top offenders as we go
	rmsd_distances_per_variant = {}
	distances_per_variant = {}
	if len(variants)!=2: panic("at the moemnt I can handle one and only one var per allele")

	# age,va = last point in progression - the age to which we want to run the simulation
	age, va = [float(number) for number in progression.split(";")[-1].split(":")]
	# expected va: simulate to  nearest rounded age
	vid = list(variants.keys())
	x, y = sim(params[vid[0]], params[vid[1]], int(round(age))+50)

	interpol = interp1d(x, y["throughput"])

	d = 0
	norm = 0
	literal_value = 0
	for agepoint in progression.split(";"):
		age, va = [float(p) for p in agepoint.split(":")]
		expected_va = interpol([age])[0]
		d += (va - expected_va)**2
		norm += 1
		literal_value += expected_va - va

	d = sqrt(d/norm) if norm>0 else 0
	literal_value = literal_value/norm  if norm>0 else 0
	parametrization_descriptor = {"rmsd": d, "sign": literal_value}
	return parametrization_descriptor

########################################
def tweak(e_prev, t_prev, dist_per_var, character):
	[cdna, prot, notes] = character
	[e,t] = [e_prev, t_prev]
	# if dist_per_var<0:
	if random() < 0:
		# we need to increase the (presumed) performance
		if 'splice' in prot:
			e = e_prev + 0.1
		else:
			if random()<0.5:
				e = e_prev + 0.1
			else:
				t = t_prev+ 0.1

		if e>1.0: e=1.0
		if t>1.0: t=1.0


	else:
		if 'splice' in prot:
			e = e_prev - 0.1
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
	for pass_number in range(1,100):
		param_descr_prev = target_function(params, case_variants, progression)
		dist_prev = param_descr_prev["rmsd"]
		print("%2d     %.2f"%(pass_number, dist_prev))

		ct = 0
		for vid, prms in params.items():
			if fixed[vid]: continue # we do not change the null variants
			ct += 1
			if not ct%50: print(f"{ct} out of {len(params)}")
			[e_prev, t_prev] = params[vid]
			[e, t] = tweak(e_prev, t_prev, param_descr_prev["sign"], character[vid])
			params[vid] = [e,t]
			param_descr = target_function(params, case_variants, progression)
			dist = param_descr["rmsd"]
			# print(" %7.2f   %7.2f   %7.2f    %7.2f   %7.2f      %7.3f   %7.3f   " %
			#       (param_descr_prev["sign"],  e_prev, t_prev, e, t, dist_prev, dist), character[vid])
			if dist<dist_prev:
				# accept
				dist_prev = dist
			else:
				# backpedal

				params[vid] = [e_prev,t_prev]

		#print("%2d     %.2f\n"%(pass_number, dist_prev))

	return params

def plot_sim_results_vs_data(x, y, progression, case_ids, title):
	age = {}
	va = {}
	for case_id in case_ids:
		age[case_id] = []
		va[case_id] = []
		for point in progression[case_id].split(";"):
			a, v = point.split(":")
			age[case_id].append(float(a))
			va[case_id].append(float(v))

	fig, ax = plt.subplots()
	plt.ylim(-0.1,1.1)
	plt.xlim(0, 30)
	ax.set_title(title)
	# ax.plot(x, y["throughput"], label='throughput')
	# ax.plot(x, y["fraction_0"], dashes=[6, 2], label='f1')
	# ax.plot(x, y["fraction_1"], dashes=[6, 2], label='f2')
	ax.plot(x, y["rpe"],  label='model')
	ax.set_xlabel("AGE")
	ax.set_ylabel("ACUITY")
	point_colors = ['red', 'blue', 'orange', 'green', 'purple', 'yellow', 'cyan', 'pink']
	for i in range(len(case_ids)):
		case_id = case_ids[i]
		ax.scatter(age[case_id], va[case_id], color=point_colors[i%len(point_colors)])
	ax.legend()

	plt.show()
	return


#########################################
def main():

	[case_variants, progression] = read_in_values()

	for case_id, vars in case_variants.items():
		parameters = {}
		character = {}
		for vid, variant_info in vars.items():
			parameters[vid] = variant_info[:2]
			character[vid]  = variant_info[2:5]

		print()
		print(case_id, f"number of vars {len(parameters)}")
		print(parameters)

		param_descr  = target_function(parameters, vars, progression[case_id], character)
		print("intial d = %.2f"% param_descr["rmsd"])
		new_variant_params = optimization_loop(parameters, vars, progression[case_id],  character)
		print(new_variant_params)
		params_vals = list(new_variant_params.values())
		print(params_vals[0],  params_vals[1])
		# simulate
		x, y = sim(params_vals[0], params_vals[1], max_age=50)
		# plot_
		# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"

		title = f"a1:  %.2f  %.2f\na2:  %.2f  %.2f" % (params_vals[0][0], params_vals[0][1], params_vals[1][0], params_vals[1][1])
		plot_sim_results_vs_data(x, y, progression, [case_id], title)


	#
	# plot(parameters, new_variant_params, character)

#########################################
if __name__ == '__main__':
	main()

