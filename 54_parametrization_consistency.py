#!/usr/bin/python3

from utils.mysql import *
from utils.abca4_mysql import  *


def variant_info(cursor, allele_id):
	qry = f"select variant_ids from alleles where id={allele_id}"
	variants = []
	for vid in hard_landing_search(cursor, qry)[0][0].strip("-").split("-"):
		qry  = "select v.cdna, v.protein,  v.gnomad_homozygotes, "
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
	qry = "select allele_id_1, allele_id_2,  onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2,  onset_age, progression] = case
		if onset_age<20: continue # the onse age not given
		vars = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue # deal with alleles with multiple variants later
			vars.extend(ret)
		if len(vars)!=2: continue
		if "Gly1961Glu" in vars[0][1] or "Gly1961Glu" in vars[1][1]: continue
		if "null" in vars[0][5] and "null" in vars[1][5]:
			print("\n onset age:", onset_age)
			print(vars)
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