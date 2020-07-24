#!/usr/bin/python3 -u


from utils.utils import  *
from utils.abca4_mysql import  *

import matplotlib.pyplot as plt

#########################################
def main():


	db, cursor = abca4_connect()

	x = []
	y = []
	fltr = "ter|splice|del|exotic|ter"

	qry  = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases "
	qry += "where onset_age>0 and (notes is null or notes not like '%caveat%')"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		[cdna1, prot1, freqs1, homozygs1, region1, cons1] = variants_from_allele(cursor, allele_id_1)
		[cdna2, prot2, freqs2, homozygs2, region2, cons2] = variants_from_allele(cursor, allele_id_2)

		# fltr = "ter|splice|del|exotic|ter"
		# if is_null(cdna1, prot1, fltr) and is_missense(prot2):
		# 	[afrm, pos, ato] = parse_protein(prot2)
		# 	x.append(pos)
		# 	y.append(onset_age)
		# elif is_null(cdna2, prot2, fltr) and is_missense(prot1):
		# 	[afrm, pos, ato] = parse_protein(prot1)
		# 	x.append(pos)
		# 	y.append(onset_age)

		fltr = "splice|del|exotic"
		if is_null(cdna1, prot1, fltr) and is_ter(cdna2, prot2):
			[afrm, pos, ato] = parse_protein(prot2)
			x.append(pos)
			y.append(onset_age)
		elif is_null(cdna2, prot2, fltr) and is_ter(cdna1, prot1):
			[afrm, pos, ato] = parse_protein(prot1)
			x.append(pos)
			y.append(onset_age)

	cursor.close()
	db.close()

	plt.stem(x,y, use_line_collection=True)

	plt.show()

#########################################
if __name__ == '__main__':
	main()
