#!/usr/bin/python3 -u


from utils.mysql import *

import xlsxwriter


def variants_from_allele(cursor, allele_id):
	qry = f"select variant_ids from alleles where id = {allele_id}"
	var_ids = [v[1:-1].split("-") for v in hard_landing_search(cursor, qry)[0]]
	for vi in var_ids:
		qry = f"select cdna, protein from variants where id in ({','.join(vi)})"
		cdna = []
		protein = []
		for [c, p] in hard_landing_search(cursor, qry):
			cdna.append(f"c.{c}" if c else "np")
			protein.append(f"p.{p}")
	return "; ".join(cdna), "; ".join(protein)


################
def write_header(worksheet, header, header_format):
	worksheet.set_row(0, 40, header_format)
	for column in range(len(header)):
		worksheet.write_string(0, column, header[column])


def table_creator(cursor, workbook, xlsx_format):
	worksheet = workbook.add_worksheet("ABCA4 gen-phen correlation")
	worksheet.set_default_row(40)
	# header
	header = ["pubmed", "reference", "patient", "allele 1: cdna",  "allele 1: protein", "allele 2: cdna",  "allele 2: protein" ]
	write_header(worksheet, header, xlsx_format["header"])
	# rows
	row = 0
	qry = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		cdna1, prot1 = variants_from_allele(cursor, allele_id_1)
		cdna2, prot2 = variants_from_allele(cursor, allele_id_2)
		row += 1
		col = 0
		for content in [publication_id, "ref", patient_xref_id, cdna1, prot1, cdna2, prot2]:
			worksheet.write(row, col, content)
			col += 1

#########################################
def main():
	db, cursor = abca4_connect()
	# Create an new Excel file and add a worksheet.
	workbook = xlsxwriter.Workbook('abca4_genphen.xlsx')
	xlsx_format = {"header":workbook.add_format({'align': 'center',  'valign': 'vcenter', 'bold': True, 'text_wrap': True}),
	               "wordwrap":workbook.add_format({'align': 'left', 'text_wrap': True}),
					"hyperlink":workbook.add_format({'align': 'center', 'color': 'blue', 'underline': 1, 'valign': 'vcenter'})}
	table_creator(cursor, workbook, xlsx_format)
	workbook.close()

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
