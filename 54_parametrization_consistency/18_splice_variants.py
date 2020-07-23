#!/usr/bin/python3

from utils.mysql import *
from utils.abca4_mysql import  *
import re

def splice_mut_pos(cdna):
	pattern = re.match("(\d+)[+\-](\d+)(\D)>(\D)", cdna)
	if not pattern: return None
	return int(pattern.group(2))

#########################################
def main():

	db, cursor = abca4_connect()
	# first position in donor
	qry = "select id, cdna from variants where cdna like '%+2%'"
	ret = error_intolerant_search(cursor, qry)
	if not ret:
		print("no ret for:\n", qry)
		exit()

	for [variant_id, cdna] in ret:
		if splice_mut_pos(cdna)!=2: continue
		for ret in hard_landing_search(cursor, f"select id from alleles where variant_ids = '-{variant_id}-'"):
			splice_allele_id = ret[0]
			qry  = "select id, allele_id_1, allele_id_1, onset_age from cases "
			qry += f"where allele_id_1={splice_allele_id} or allele_id_2={splice_allele_id}"
			print(variant_id, cdna)
			for [case_id, allele_id_1, allele_id_2, onset_age] in hard_landing_search(cursor, qry):
				the_other_allele_id = allele_id_1 if allele_id_2==splice_allele_id else allele_id_2

				qry = f"select variant_ids from alleles where id={the_other_allele_id}"
				variant_ids = hard_landing_search(cursor, qry)[0][0].strip("-").split("-")
				if len(variant_ids)>1: continue

				qry = f"select notes from parametrization where variant_id={variant_ids[0]}"
				notes = hard_landing_search(cursor, qry)[0][0]
				if not "fs" in notes and not "stop" in notes: continue
				print("\t", case_id, onset_age, notes)

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