#!/usr/bin/python3

from sys import argv
from utils.abca4_mysql import  *
from utils.utils import unpack_progression
from utils.simulation import *
import matplotlib.pyplot as plt



def variant_params(cursor, variant_ids):
	params = {}
	te = None
	for vid in variant_ids:
		for parm_table in ['parametrization_adjusted', 'parametrization_literature', 'parametrization']:
			qry  = "select expression_folding_membrane_incorporation, transport_efficiency  "
			qry += f"from {parm_table} where variant_id={vid}"
			ret = error_intolerant_search(cursor, qry)
			if not ret: continue
			if len(ret)>2:
				print(f"multiple entries fo variant {vid} (?!)")
				exit()
			print(f"parms for variant id {vid} found in {parm_table}")
			params[vid] = ret[0]
			break
	return params

def read_baseline(cursor, case_id):
	baseline = 0.1
	qry = f"select baseline from case_parametrization where case_id={case_id}"
	ret =  error_intolerant_search(cursor, qry)
	if ret: baseline = ret[0][0]
	return baseline

def read_params(cursor, case_id, variant_ids):
	return  variant_params(cursor, variant_ids), read_baseline(cursor, case_id)

####################################
def read_progression(cursor, case_id):
	# fitting on our own cases on purpose 52 is our set of data
	qry  = "select allele_id_1, allele_id_2,  progression, erg_r_rod from cases "
	qry += f"where id = {case_id}"

	[allele_id_1, allele_id_2, progression, erg_r_rod] = hard_landing_search(cursor, qry)[0]
	variants = []
	for ai in [allele_id_1, allele_id_2]:
		qry = f"select variant_ids from alleles where id = {ai}"
		ret = hard_landing_search(cursor, qry)[0][0].strip('-').split('-')
		variants.extend(ret)

	if len(variants)>2: return ['error','multiple variants', variants] # deal with alleles with multiple variants later

	return [variants, progression, erg_r_rod]




#########################################
def report(parameters, rpe_baseline, new_variant_params, new_rpe_baseline, var_ids, case_id):
	print("\n========================")
	print(f"case id {case_id}")

	print( f"{var_ids[0]}:  %.2f  %.2f\n{var_ids[1]}:  %.2f  %.2f\nrpe_baseline:   %.2f" % \
	       (parameters[var_ids[0]][0], parameters[var_ids[0]][1],
	        parameters[var_ids[1]][0], parameters[var_ids[1]][1],  rpe_baseline[case_id]))
	print("--------------------------")

	print( f"{var_ids[0]}:  %.2f  %.2f\n{var_ids[1]}:  %.2f  %.2f\nrpe_baseline:   %.2f" % \
	      (new_variant_params[var_ids[0]][0], new_variant_params[var_ids[0]][1],
	       new_variant_params[var_ids[1]][0], new_variant_params[var_ids[1]][1],  new_rpe_baseline[case_id]))
	print()

def plot_sim_results_vs_data(age_va, va, age_erg, erg, varid1, varid2, params1, params2, rpe_baseline):

	# simulate
	x, y = sim(params1, params2, rpe_baseline, max_age=50)
	# plot_
	# title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
	title = f"{varid1}:  %.2f  %.2f\n{varid2}:  %.2f  %.2f" % (params1[0], params1[1], params2[0], params2[1])

	fig, axs = plt.subplots()
	axs.set_ylim(-0.1,1.1)
	axs.set_title(title)
	axs.plot(x, y["rpe"])
	axs.scatter(age_va, va, color='red')
	axs.scatter(age_erg, erg, color='green')

	plt.show()
	return

def visualization(case_id, varids, params, baseline, va_progression, erg_progression):

	# varids = list(case_variants[case_id].keys())
	# report(old_params, old_rpe_baseline, params, rpe_baseline, varids, case_id)
	age_va, va = unpack_progression(va_progression)
	age_erg, erg = unpack_progression(erg_progression)
	erg_norm = 400.0
	erg = [float("%.2f"%(e/erg_norm)) for e in erg]

	plot_sim_results_vs_data(age_va, va, age_erg, erg, varids[0], varids[1],
	                         params[varids[0]],  params[varids[1]], baseline)

	pass


#########################################
def process_case(cursor, case_id):

	[case_variants, progression, erg] = read_progression(cursor, case_id)
	if 'error' in case_variants:
		print(case_variants, progression, erg)
		return
	[params, baseline] = read_params(cursor, case_id, case_variants)
	print("case variants:", case_variants)
	print("params:", params)
	print("baseline:", baseline)
	print("progression:", progression)
	print("erg:", erg)

	visualization(case_id, case_variants, params, baseline, progression, erg)

#########################################
def main():
	db, cursor = abca4_connect()
	if len(argv) >= 2:
		case_ids = [argv[1]]
	else:
		# print(f"usage {argv[0]} <case_id>")
		# exit()
		qry = "select id from cases where erg_r_rod is not null"
		case_ids = [r[0] for r in hard_landing_search(cursor, qry)]
	for case_id in case_ids:
		process_case(cursor, case_id)
	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

