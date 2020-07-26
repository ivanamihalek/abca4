#!/usr/bin/python3

'''
Mark as supicious cases where:
1) one allele has single synonymous variant (which is not at the splice position)
2) cases with two null alleles with onset after 15 yrs
3) one allele has a single deep intronic variant with no experimental support
4) Gly1961Glu
5) more than 10 homozygotes in gnomAD
'''

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
	qry = "select id, allele_id_1, allele_id_2,  onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[case_id, allele_id_1, allele_id_2,  onset_age, progression] = case

		vars = []
		for ai in [allele_id_1, allele_id_2]:
			ret = variant_info(cursor, ai)
			if len(ret)>2: continue
			vars.extend(ret)
		if len(vars)!=2: continue
		notes = []
		for i in range(2):
			mutation =  vars[i][1]
			if "Gly1961Glu" in mutation:
				notes.append(f"caveat allele {i+1}: Gly1961Glu")
				continue
			homozygotes = vars[i][2]
			if homozygotes and homozygotes>10:
				notes.append(f"caveat allele {i+1}: >10 homozygotes")
				continue

			variant_notes =  vars[i][5]
			if "suspicious" in variant_notes: # these are, so far deep intronic and
				notes.append(f"caveat allele {i+1}: {variant_notes.replace('suspicious: ', '')}")
		if "Ter" in vars[0][1] and "Ter" in vars[1][1] and onset_age>15:
			notes.append("caveat: late onset for double null")
		if notes:
			notestring = "; ".join(notes)
			qry = f"update cases set notes='{notestring}' where id={case_id}"
			print(qry)
			error_intolerant_search(cursor, qry)
	db.close()


#########################################
if __name__ == '__main__':
	main()
