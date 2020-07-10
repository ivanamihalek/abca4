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
def column_string(idx):
	char = chr(ord('A')+idx)
	return f"{char}:{char}"

def set_column_widths(worksheet, header, wwrap_format):
	# we'll put image in the second column - not sure what are the units here
	# here:  https://stackoverflow.com/questions/47345811/excel-cell-default-measure-unit
	# says that One unit of column width is equal to the width of one character in the Normal style (?)

	#idx = header.index("protein position")
	#worksheet.set_column(column_string(idx), len("position"))

	for title in header:
		idx = header.index(title)
		worksheet.set_column(column_string(idx), len(title)+5)

	for title in ["reference"]:
		idx = header.index(title)
		worksheet.set_column(column_string(idx), 30, wwrap_format)

	for title in [ "allele 1: protein",  "allele 2: protein"]:
		idx = header.index(title)
		worksheet.set_column(column_string(idx), 2*len(title), wwrap_format)


def write_header(worksheet, header, header_format):
	worksheet.set_row(0, 40, header_format)
	for column in range(len(header)):
		worksheet.write_string(0, column, header[column])


def table_creator(cursor, workbook, xlsx_format):
	worksheet = workbook.add_worksheet("ABCA4 gen-phen correlation")
	worksheet.set_default_row(25)

	# header
	header = ["pubmed", "reference", "patient", "allele 1: cdna",  "allele 1: protein", "allele 2: cdna",  "allele 2: protein" ]
	set_column_widths(worksheet, header, xlsx_format["wordwrap"])
	write_header(worksheet, header, xlsx_format["header"])
	# rows
	row = 0
	qry = "select allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression from cases"
	for case in hard_landing_search(cursor, qry):
		[allele_id_1, allele_id_2, publication_id, patient_xref_id, onset_age, progression] = case
		cdna1, prot1 = variants_from_allele(cursor, allele_id_1)
		cdna2, prot2 = variants_from_allele(cursor, allele_id_2)
		row += 1
		column = 0

		qry = f"select pubmed, pubmedcentral, reference  from publications where id = {publication_id}"
		[pubmed_id, pmc_id, reference] = hard_landing_search(cursor,qry)[0]
		pubmed_hyperlink = "http://pubmed.ncbi.nlm.nih.gov/%s" % pubmed_id
		worksheet.write_url(row, column, pubmed_hyperlink, string=str(pubmed_id))
		column += 1

		pmc_hyperlink = "https://www.ncbi.nlm.nih.gov/pmc/articles/%s" % pmc_id
		worksheet.write_url(row, column, pmc_hyperlink, string=reference)
		column += 1

		for content in [patient_xref_id, cdna1, prot1, cdna2, prot2]:
			worksheet.write(row, column, content)
			column += 1

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
