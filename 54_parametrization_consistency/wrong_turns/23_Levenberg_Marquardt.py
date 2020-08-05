#!/usr/bin/python3
import copy
from math import sqrt, exp
from random import random

from utils.abca4_mysql import  *

from utils.simulation import progression_sim
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import least_squares

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

def sim(var1_param, var2_param, rpe_baseline, max_age):
	# the values for one case/patient
	# alpha is the abilty to fold and incorporate
	# "fraction" refers to the fraction of the wild-type capability
	expression_fraction = [var1_param[0], var2_param[0]]
	transport_fraction = [var1_param[1], var2_param[1]]
	x, y  = progression_sim(max_age, expression_fraction, transport_fraction, rpe_baseline=rpe_baseline, verbose=False)
	return x, y


#########################################
# turn  the target into something that scipy understands
def target_function(parameters, age, va):
	params_variant_1 = parameters[:2]
	params_variant_2 = parameters[2:4]
	rpe_baseline = parameters[4]

	max_age = int(age[-1]) + 2
	x, y = sim(params_variant_1, params_variant_2, rpe_baseline, max_age)
	lower_intrpl_point = max(int(age[0])-1, 0)
	upper_intrpl_point = max_age
	intrpl_range = x[lower_intrpl_point:upper_intrpl_point]
	intrpl_values =  y["rpe"][lower_intrpl_point:upper_intrpl_point]
	interpol = interp1d(intrpl_range,intrpl_values)
	residuals = []

	return residuals


def plot_sim_results_vs_data(age, va, params, new_params):

	# simulate
	x, y = sim(params[0:2], params[2:4], params[4], max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	title = f"a1:  %.2f  %.2f\na2:  %.2f  %.2f" % (params[0], params[1], params[2], params[3])

	# simulate
	new_x, new_y = sim(new_params[0:2], new_params[2:4], new_params[4], max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	new_title = f"a1:  %.2f  %.2f\na2:  %.2f  %.2f" % (new_params[0], new_params[1], new_params[2], new_params[3])

	rows = 1
	cols = 2
	fig, axs = plt.subplots(rows, cols)
	axs[0].set_title(title)
	axs[0].plot(x, y["rpe"])
	axs[0].scatter(age, va, color='red')

	axs[1].set_title(new_title)
	axs[1].plot(new_x, new_y["rpe"])
	axs[1].scatter(age, va, color='red')

	plt.show()
	return


#########################################
def unpack_progression(progression):
	age = []
	va = []
	for agepoint in progression.split(";"):
		[a, v]  = [float(number) for number in agepoint.split(":")]
		age.append(a)
		va.append(v)
	return age, va

def eyeball_bounds(parameters):
	lower = []
	upper = []
	for p in parameters[:4]:
		lower.append(max(p-0.25, 0))
		upper.append(min(p+0.25, 1))
	lower.append(0)
	upper.append(0.2)

	return [lower, upper]

##########
def main():

	print("this worx wors then my MC walk")
	exit()

	[case_variants, progression] = read_in_values()
	rpe_baseline = 0.1


	for case_id, variants in case_variants.items():

		print()
		print(case_id, f"number of vars {len(variants)}")

		parameters = []
		character  = {}
		for vid, variant_info in variants.items():
			parameters.extend(variant_info[:2])
			character[vid] = variant_info[2:5]
		parameters.append(rpe_baseline)
		parameter_bounds = eyeball_bounds(parameters)
		print(character)
		print("initial params", parameters)

		age, va = unpack_progression(progression[case_id])
		print("initial residuals", ["%.2f"%res for res in target_function(parameters, age, va)])

		lsq_return = least_squares(target_function, parameters, bounds=parameter_bounds,  args=(age, va), method='trf')
		# lsq_return is pretty much a dictionary;
		# intersting keys are message, fun (which actually returns new residuals)
		# and x which returns optimized parameter values
		new_params = lsq_return.get("x").tolist()
		new_residuals =  lsq_return.get("fun").tolist()
		print("new params",["%.3f"%p for p in  new_params])
		print("new residuals         ", ["%.2f"%res for res in target_function(new_params, age, va)])
		print("new residuals reported", ["%.2f"%res for res in new_residuals])

		plot_sim_results_vs_data(age, va, parameters, new_params)

		exit()

	#
	# plot(parameters, new_variant_params, character)

#########################################
if __name__ == '__main__':
	main()

