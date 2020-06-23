#!/usr/bin/python3

import os
import re
import subprocess

from utils.mysql import *
from utils.utils import single_letter_code
from sys import argv
from Bio.Seq  import Seq
from Bio.Alphabet import generic_dna, generic_protein

def get_canonical_exon_bdries(cursor, ensembl_gene_id):

	qry = f"select canonical_transcript_id, seq_region_id, seq_region_strand from gene where stable_id='{ensembl_gene_id}'"
	[canonical_transcript_id, seq_region_id, seq_region_strand] =  hard_landing_search(cursor, qry)[0]
	print(canonical_transcript_id, seq_region_id, seq_region_strand)
	neg_strand = seq_region_strand<0

	# translation start and end
	qry  = "select seq_start, start_exon_id,  seq_end, end_exon_id "
	qry += f" from translation where transcript_id = {canonical_transcript_id} "
	[transl_start, start_exon_id, transl_end, end_exon_id] = hard_landing_search(cursor, qry)[0]

	# exons
	qry = f"select exon_id from exon_transcript  where transcript_id = {canonical_transcript_id}"
	exon_ids = [row[0] for row in hard_landing_search(cursor, qry)] # dict giving exon_id as function of rank

	coordinates = []
	reading = False
	for exon_id in exon_ids:
		qry = f"select  seq_region_start, seq_region_end from exon where exon_id = {exon_id}"
		[seq_region_start, seq_region_end] = hard_landing_search(cursor, qry)[0]
		if exon_id==start_exon_id:
			reading = True
			seq_region_end -= transl_end
			print("start", seq_region_start, transl_start)
		if exon_id==end_exon_id:
			reading = False
			seq_region_start += transl_end

			print("end", seq_region_end, transl_end)
		if not reading: continue
		coordinates.append([seq_region_start, seq_region_end])


	total_length = 0
	bdries = [0]
	coordinates.sort(key=lambda x: x[0])
	for [seq_region_start, seq_region_end] in coordinates:
		length = seq_region_end-seq_region_start+1
		print( seq_region_start, seq_region_end, "   %5d"%length)
		total_length += length
		bdries.append(total_length)


	print(total_length, total_length%3, total_length/3)
	# print(bdries)
	return seq_region_id, seq_region_strand, coordinates, bdries


#########################################
def get_cdna(cursor, blastdbcmd, cdna_fasta, seq_region_id, strand,  raw_coordinates):
	print(strand)
	qry =  f"select name from seq_region where seq_region_id={seq_region_id}"
	seq_region_name = hard_landing_search(cursor, qry)[0][0]

	seq = ""
	for [start, end] in raw_coordinates:
		tmpfile = "tmp.fasta"
		if os.path.exists(tmpfile): os.remove(tmpfile)
		cmd  = f"{blastdbcmd} -db {cdna_fasta} -dbtype nucl -entry {seq_region_name} "
		cmd += f"-range {start}-{end} -out {tmpfile} -outfmt %s"
		subprocess.call(["bash", "-c", cmd])
		if not os.path.exists(tmpfile):
			print(f"{tmpfile} not produced")
			exit()
		with open(tmpfile) as inf:
			seq += inf.read().replace("\n", "")
	print(len(seq))
	biopython_dna = Seq(seq, generic_dna).reverse_complement()
	print(biopython_dna)
	print("============")
	print(biopython_dna.translate())
	exit()

########################################
def main():

	blastdbcmd = "/usr/bin/blastdbcmd"
	cdna_fasta = "/storage/databases/ensembl-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.fa"
	for dep in [blastdbcmd, cdna_fasta]:
		if not os.path.exists(dep):
			print(f"{dep} not found")
			exit()

	db,cursor = abca4_connect()

	# ensmebl id from gene name
	qry = "select ensembl_gene_id from identifier_maps.hgnc where approved_symbol='ABCA4'"
	ensembl_gene_id = hard_landing_search(cursor, qry)[0][0]
	print(ensembl_gene_id)

	switch_to_db(cursor, "homo_sapiens_core_97_38")

	# coordinates of canonical exons
	# we will just extract positions of exon boundaries on cdna
	[seq_region_id, strand, raw_coordinates, canonical_exon_bdries] = get_canonical_exon_bdries(cursor, ensembl_gene_id)

	# canonical cdna
	cdna = get_cdna(cursor, blastdbcmd, cdna_fasta, seq_region_id, strand, raw_coordinates)

	# canonical translation

	# make sure they match

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

