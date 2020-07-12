#!/usr/bin/python3 -u


from utils.abca4_mysql import  *
import re
import matplotlib.pyplot as plt

#########################################
def main():

	db, cursor = abca4_connect()

	qry = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		print()
		print("================================")
		print(f"{publication_id}    {patient_xref_id}               {onset_age}     {progression}")
		for aid in [allele_id_1, allele_id_2]:
			print(f"\t allele  {aid}")
			qry = f"select variant_ids from alleles where id = {aid}"
			var_ids = hard_landing_search(cursor, qry)[0][0][1:-1].split("-")
			for vid in var_ids:
				qry = f"select cdna, protein, gnomad_freq, gnomad_homozygotes from variants where id={vid}"
				[cdna, protein, gnomad_freq, gnomad_homozygote] =  hard_landing_search(cursor, qry)[0]
				print(f"\t\t   variant  {vid}   {cdna}   {protein}  {gnomad_freq}   {gnomad_homozygote} ")


	cursor.close()
	db.close()
#########################################
if __name__ == '__main__':
	main()
