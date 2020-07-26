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
		qry += "from variants v, parametrization p "
		qry += f"where p.variant_id=v.id and v.id={vid}"
		ret = hard_landing_search(cursor, qry)
		if len(ret)>2:
			print(f"multiple entries fo variant {vid} (?!)")
			exit()
		variants.append(ret[0])
	return variants


def read_in_values():
	db, cursor = abca4_connect()
	qry  = "select id, allele_id_1, allele_id_2,  progression from cases "
	qry += "where (notes is null or notes not like '%caveat%') "
	# we don't want to parametrize on our own cases, bu they are already filtered out because they do not contain the onset age
	# (the onset age is set to -1)
	prog = {}
	case_variants = {}
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2, progression ] = case
		#if len(progression.split(";"))<3: continue
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

def plot_sim_results_vs_data(x, y,  progression, title):
	age = []
	va = []
	for point in progression.split(";"):
		a, v = point.split(":")
		age.append(float(a))
		va.append(float(v))

	fig, ax = plt.subplots()
	plt.ylim(0,1.1)
	#plt.xlim(0, int(max(age)+2))
	ax.set_title(title)
	ax.plot(x, y["throughput"], label='throughput')
	ax.plot(x, y["fraction_0"], dashes=[6, 2], label='f1')
	ax.plot(x, y["fraction_1"], dashes=[6, 2], label='f2')
	ax.plot(x, y["rpe"],  label='RPE')
	ax.scatter(age, va, color='red', label='patient')
	ax.legend()

	plt.show()
	return


def main ():
	# read
	[case_variants, progression] = read_in_values()
	for case_id, variants in case_variants.items():
		print(case_id, variants)
		[e1, t1, notes1, cdna1, prot1, vid1]  = variants[0]
		[e2, t2, notes2, cdna2, prot2, vid2]  = variants[1]
		# simulate
		x, y = sim([e1, e2], [t1, t2])
		# plot_
		title = f"a1, {vid1}: {cdna1} {prot1} {e1} {t1} \na2, {vid2}: {cdna2} {prot2} {e2} {t2}"
		plot_sim_results_vs_data(x, y, progression[case_id], title)

####################
if __name__ == "__main__":
	main()
