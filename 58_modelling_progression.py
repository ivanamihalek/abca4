#! /usr/bin/python3

from utils.abca4_mysql import  *

from utils.simulation import progression_sim
import matplotlib.pyplot as plt

def variant_info(cursor, allele_id):
	qry = f"select variant_ids from alleles where id={allele_id}"
	variants = []
	for vid in hard_landing_search(cursor, qry)[0][0].strip("-").split("-"):
		qry  = "select  p.expression_folding_membrane_incorporation, p.transport_efficiency, p.notes, "
		qry += "v.cdna, v.protein, v.id "
		qry += "from variants v, parametrization_adjusted p "
		qry += f"where p.variant_id=v.id and v.id={vid}"
		ret = error_intolerant_search(cursor, qry)
		if not ret: continue # we are not dealing with this varaint for some reason
		if len(ret)>2:
			print(f"multiple entries fo variant {vid} (?!)")
			exit()
		variants.append(ret[0])
	return variants


def read_in_values():
	db, cursor = abca4_connect()
	qry  = "select id, allele_id_1, allele_id_2,  progression from cases "
	qry += "where (notes is null or notes not like '%caveat%')  and progression is not null and progression!=''"
	# we don't want to parametrize on our own cases, bu they are already filtered out because they do not contain the onset age
	# (the onset age is set to -1)
	prog = {}
	case_variants = {}
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2, progression ] = case
		if len(progression.split(";"))<3: continue
		variants = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			variants.extend(ret)
		if len(variants)!=2: continue

		case_variants[case_id] = variants
		prog[case_id] = progression

	cursor.close()
	db.close()

	return [case_variants, prog]


def sim(expression_fraction, transport_fraction):
	# the values for one case/patient
	# alpha is the abilty to fold and incorporate
	# "fraction" refers to the fraction of the wild-type capability
	max_age = 80
	x, y  = progression_sim(max_age, expression_fraction, transport_fraction, verbose=False)
	return x, y

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


def main ():
	# read
	[case_variants, progression] = read_in_values()

	common_genotype = {}
	for case_id, variants in case_variants.items():
		#print(case_id, variants)
		[e1, t1, notes1, cdna1, prot1, vid1]  = variants[0]
		[e2, t2, notes2, cdna2, prot2, vid2]  = variants[1]
		genotype = "_".join([str(v) for v in sorted([vid1, vid2])])
		if genotype not in common_genotype: common_genotype[genotype] = []
		common_genotype[genotype].append(case_id)

		# simulate
		x, y = sim([e1, e2], [t1, t2])
		# plot_
		title = f"a1: {cdna1} {prot1} {e1} {t1} \na2: {cdna2} {prot2} {e2} {t2}"
		plot_sim_results_vs_data(x, y, progression, [case_id], title)

	# for genotype, case_ids in common_genotype.items():
	# 	if len(case_ids)<2: continue
	# 	print(genotype, case_ids)
	# 	case_id = case_ids[0]
	# 	variants = case_variants[case_id] # they should all be the same, right
	# 	[e1, t1, notes1, cdna1, prot1, vid1]  = variants[0]
	# 	[e2, t2, notes2, cdna2, prot2, vid2]  = variants[1]
	# 	x, y = sim([e1, e2], [t1, t2])
	# 	# plot_
	# 	title = f"a1, {vid1}: {cdna1} {prot1} {e1} {t1} \na2, {vid2}: {cdna2} {prot2} {e2} {t2}"
	# 	plot_sim_results_vs_data(x, y, progression, case_ids, title)

####################
if __name__ == "__main__":
	main()
