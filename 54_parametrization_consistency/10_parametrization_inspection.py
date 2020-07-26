#!/usr/bin/python3

from utils.mysql import *
from utils.abca4_mysql import  *
from math import sqrt
from random import random
from statistics import mean, stdev

import matplotlib.pyplot as plt

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
	qry += "where onset_age>0 and (notes is null or notes not like '%caveat%') "
	# we don't want to parametrize on our own cases, bu they are already filtered out because they do not contain the onset age
	# (the onset age is set to -1)
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

	cursor.close()
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


##############
def distros(case_variants, age, variant_parameters):

	db, cursor = abca4_connect()

	age_list = {}
	cases_by_label = {}
	for case_id, vars in case_variants.items():
		label = " | ".join(sorted([variant_parameters[v] for v in vars]))
		if not label in age_list:
			age_list[label] = []
			cases_by_label[label] = []
		age_list[label].append(age[case_id])
		cases_by_label[label].append(case_id)


	labels_sorted = sorted(age_list.keys())
	for label in labels_sorted:
		agelist  = age_list[label]
		m = mean(agelist)
		s = stdev(agelist) if len(agelist)>1 else 0
		print(label, "       %2.0f     %2.0f"%(m, s),  "      ", sorted(agelist))
		#if label == "1_4 | 1_4":
		if s>0:
			outlayers = set(filter(lambda x: abs(x-m)/s>3 or x<6, agelist))
			for onset_age in outlayers:
				case_ids = list(filter(lambda case_id: age[case_id]==onset_age, cases_by_label[label]))
				for cid in case_ids:
					print("\t\t ", onset_age, cid)
					for vid in case_variants[cid]:
						qry  = "select cdna, protein, gnomad_homozygotes, conserved_in_ortho_verts, conserved_in_verts_insects "
						qry += f"from variants where id={vid}"
						print("\t\t\t", hard_landing_search(cursor, qry)[0])
		print()
		# exit()

	print(f"number of labels {len(age_list)}")
	cursor.close()
	db.close()

	return age_list


def label2name(label):
	name = ""
	vct = 0
	for v in label.split(" | "):
		vct +=1
		[e, t] = ["%3d"%(25*int(i)) for i in  v.split("_")]
		if e== "%3d"%0:
			name  += f"var{vct}: null\n"
		else:
			name  += f"var{vct}: e={e}%  t={t}%\n"
	return name

def plot_all_distros(age_list):

	number_of_bins = 75
	bin_size = 1
	bins = []
	for b in range(number_of_bins):
		bins.append(b*bin_size+0.1)

	rows = 5
	cols = 4
	fig, axs = plt.subplots(rows, cols)
	fig.text(0.5, 0.04, 'onset age', ha='center')
	fig.text(0.04, 0.5, 'number of patients', va='center', rotation='vertical')

	labels_sorted = sorted(age_list.keys())
	n = -1
	for label in labels_sorted:
		agelist  = age_list[label]
		if len(agelist)<5: continue
		n += 1
		if n>=rows*cols: break
		name = label2name(label)
		axs.flat[n].hist(agelist, bins)
		#axs.flat[n].legend([f"{name} ({len(agelist)})"])
		axs.flat[n].text(0.6, 0.75, f"{name}", horizontalalignment='center',
		                 verticalalignment='center',  transform = axs.flat[n].transAxes)
		axs.flat[n].axvline(x=5, color="red")
		axs.flat[n].axvline(x=10, color="red")

	plt.show()


def plot_distros_with_null(age_list):

	number_of_bins = 75
	bin_size = 1
	bins = []
	for b in range(number_of_bins):
		bins.append(b*bin_size+0.1)

	rows = 3
	cols = 3
	fig, axs = plt.subplots(rows, cols)
	fig.text(0.5, 0.04, 'onset age', ha='center')
	fig.text(0.04, 0.5, 'number of patients', va='center', rotation='vertical')

	n = 0
	label = f"0_4 | 0_4"
	name = label
	agelist = age_list[label]
	axs.flat[n].hist(agelist, bins)
	axs.flat[n].legend([f"{name} ({len(agelist)})"])
	axs.flat[n].axvline(x=5, color="red")
	axs.flat[n].axvline(x=10, color="red")

	other_alleles =  [       "1_4", "2_4",
	                  "4_1", "2_2", "3_4",
	                  "4_2", "4_3", "3_3"]

	for other_allele in other_alleles:
		n += 1
		if n>=rows*cols: break
		label = f"0_4 | {other_allele}"
		print(label)
		agelist = age_list.get(label, [])
		name = label
		axs.flat[n].ylim(0,10)
		axs.flat[n].hist(agelist, bins)
		#axs.flat[n].legend([f"{name} ({len(agelist)})"])
		axs.flat[n].text(number_of_bins/2, 7, f"{name}")
		axs.flat[n].axvline(x=5, color="red")
		axs.flat[n].axvline(x=10, color="red")

	plt.show()


#########################################
def main():

	[case_variants, age, variant_parameters] = read_in_values()
	age_list = distros(case_variants, age, variant_parameters)
	# plot_distros_with_null(age_list)
	plot_all_distros(age_list)

#########################################
if __name__ == '__main__':
	main()

