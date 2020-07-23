#!/usr/bin/python3

from utils.mysql import *
from utils.abca4_mysql import  *


def variant_info(cursor, allele_id):
	qry = f"select variant_ids from alleles where id={allele_id}"
	variants = []
	for vid in hard_landing_search(cursor, qry)[0][0].strip("-").split("-"):
		qry  = "select v.id, v.cdna, v.protein,  v.gnomad_homozygotes, "
		qry += "p.expression_folding_membrane_incorporation, p.transport_efficiency, p.notes "
		qry += "from variants v, parametrization p "
		qry += f"where p.variant_id=v.id and v.id={vid}"
		ret = hard_landing_search(cursor, qry)
		if len(ret)>2:
			print(f"multiple entries fo variant {vid} (?!)")
			exit()
		variants.append(ret[0])
	return variants

#########################################
# if it appears that a mutant should be null, consdire distribution in other cases
# do other cases have an accompanying null?
def putative_null(cursor, query_case_id, variant_id):
	qry  = f"select id from alleles where variant_ids = '-{variant_id}-'"
	aids = [str(ret[0]) for ret in hard_landing_search(cursor, qry)]
	aidstring = ",".join(aids)

	qry  = "select id, allele_id_1, allele_id_2,  onset_age, progression from cases "
	qry += f"where onset_age>0  and (allele_id_1 in ({aidstring}) or allele_id_2 in ({aidstring})) "
	qry += f"and id!={query_case_id} and (notes is null or notes not like '%caveat%')"
	ret = error_intolerant_search(cursor, qry)
	if not ret: return
	for case in sorted(ret, key=lambda x: x[3]): # sort by age of onset
		[case_id, allele_id_1, allele_id_2,  onset_age, progression] = case
		variants = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			variants.extend(ret)
		print(f"\nonset age: {onset_age}  case_id:{case_id}")
		print(variants)

	# exit()

#########################################
def main():

	db, cursor = abca4_connect()
	qry  = "select id, allele_id_1, allele_id_2,  onset_age, progression from cases "
	qry += "where onset_age>0 and (notes is null or notes not like '%caveat%')"
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2,  onset_age, progression] = case
		if onset_age>9: continue
		variants = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			variants.extend(ret)
		if len(variants)!=2: continue
		if "null" in variants[0][6] and not "null" in variants[1][6]:
			print("\n\n====================================================")
			print(f"\nonset age:{onset_age}  case_id:{case_id}")
			print(variants)
			putative_null(cursor, case_id, variants[1][0])
		if "null" in variants[1][6] and not "null" in variants[0][6]:
			print("\n\n====================================================")
			print(f"\nonset age:{onset_age}  case_id:{case_id}")
			print(variants)
			putative_null(cursor,  case_id, variants[0][0])
	db.close()


#########################################
if __name__ == '__main__':
	main()


'''
is parametrization consistent with the pehontype produced by variants (are come variants effectively null?)
then run simulation for the initial parametrization and see how far off we are from the reported acuity of vision
does my cire sim represent a single day?
read in our patients to the database

'''