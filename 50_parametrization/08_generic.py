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
def main():

	db, cursor = abca4_connect()

	#case_ids = "139,420,432,447,449,453,458,463,482,483,584"
	#case_ids = "2,13,15,24,30,44,48,56,87,110,129,130,140,201,229,231,251,252,309,311,321,343,344,345,349,357,393,416,448,636,637"

	qry  = "select id, allele_id_1, allele_id_2,  onset_age, progression from cases "
	qry += f"where id in ({case_ids})"

	cases =  hard_landing_search(cursor, qry)
	for case in sorted(cases, key=lambda c: c[3]):
		[case_id, allele_id_1, allele_id_2,  onset_age, progression] = case
		variants = []
		print()
		print(f" {case_id}       {onset_age}     {progression}")
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			print(ret)
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