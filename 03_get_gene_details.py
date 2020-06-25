#!/usr/bin/python3

# somebody somehwere wrote a buggy program and its output became a norm:
# sometimes, but no systematically, for a gene coded in reverse strand
# when mutations are within the exons, their position is numbere from left to right

# c.3607G>T means that the position 3607 on the coding DNA (counted from the left) is mutated to T
# however, c.3523-2A>G means that at position  3523 counted from *right* there is mutation in the intron,
# 2 positions to the right(?) whichever imbecile came up with that

import os
import subprocess

from utils.mysql import *
from Bio.Seq  import Seq
from Bio.Alphabet import generic_dna, generic_protein
import pprint

class Exon:
	def __init__ (self, exon_id, seq_region_start, seq_region_end):
		self.id = exon_id
		self.start = seq_region_start
		self.end = seq_region_end
		self.is_first = False
		self.is_last = False
		self.donor_splice = None
		self.acceptor_splice = None

def get_canonical_exon_bdries(cursor, ensembl_gene_id):

	qry = f"select canonical_transcript_id, seq_region_id, seq_region_strand from gene where stable_id='{ensembl_gene_id}'"
	[canonical_transcript_id, seq_region_id, seq_region_strand] =  hard_landing_search(cursor, qry)[0]
	# print(canonical_transcript_id, seq_region_id, seq_region_strand)
	neg_strand = seq_region_strand<0
	# translation start and end - they are never 0 in the original ensembl entries, so we take the offset is 1
	qry  = "select seq_start, start_exon_id,  seq_end, end_exon_id "
	qry += f" from translation where transcript_id = {canonical_transcript_id} "
	[transl_start, start_exon_id, transl_end, end_exon_id] = hard_landing_search(cursor, qry)[0]

	# exons
	qry = f"select exon_id from exon_transcript  where transcript_id = {canonical_transcript_id}"
	exon_ids = [row[0] for row in hard_landing_search(cursor, qry)] # dict giving exon_id as function of rank


	exons = []
	for exon_id in exon_ids:
		qry = f"select  seq_region_start, seq_region_end from exon where exon_id = {exon_id}"
		[seq_region_start, seq_region_end] = hard_landing_search(cursor, qry)[0]
		new_exon = Exon(exon_id, seq_region_start, seq_region_end)
		if neg_strand:
			if exon_id==end_exon_id:
				new_exon.is_first = True
				new_exon.start = new_exon.end - (transl_end - 1)
			if exon_id==start_exon_id:
				new_exon.is_last = True
				new_exon.end -= transl_start-1
		else:
			if exon_id==end_exon_id:
				new_exon.is_last = True
				new_exon.end = new_exon.start + transl_end - 1
			if exon_id==start_exon_id:
				new_exon.is_first = True
				new_exon.start += transl_start-1

		exons.append(new_exon)

	exons.sort(key=lambda x: x.start)

	total_length = 0
	reading = True
	for exon in exons:
		if exon.is_first: reading=True
		if reading:
			length = exon.end - exon.start + 1
			total_length += length
		if exon.is_last: reading=False

	if total_length%3 != 0:
		print("total length not divisible by 3")
		exit()
	# print(total_length, total_length%3, total_length/3)
	# print()
	return seq_region_id, seq_region_strand, exons


#########################################
def get_cdna(cursor, blastdbcmd, cdna_fasta, seq_region_id, strand,  exons):

	splice_length = 12 # how many places should I put there? is this meaningful

	qry =  f"select name from seq_region where seq_region_id={seq_region_id}"
	seq_region_name = hard_landing_search(cursor, qry)[0][0]
	reverse = strand<0

	seq = ""
	for exon in exons:
		tmpfile = "tmp.fasta"
		if os.path.exists(tmpfile): os.remove(tmpfile)
		cmd  = f"{blastdbcmd} -db {cdna_fasta} -dbtype nucl -entry {seq_region_name} "
		cmd += f"-range {exon.start-splice_length}-{exon.end+splice_length} -out {tmpfile} -outfmt %s"
		subprocess.call(["bash", "-c", cmd])
		if not os.path.exists(tmpfile):
			print(f"{tmpfile} not produced")
			exit()
		with open(tmpfile) as inf:
			inseq = inf.read().replace("\n", "")
		seq += inseq[splice_length:-splice_length]
		if reverse:
			exon.donor_splice    = str(Seq(inseq[:splice_length]).reverse_complement())
			exon.acceptor_splice = str(Seq(inseq[-splice_length:]).reverse_complement())
		else:
			# if the exon is not the first, technically this is not the splice region
			exon.acceptor_splice = inseq[:splice_length]
			exon.donor_splice    = inseq[-splice_length:]
		# I do not know wtf this is, but it correposnds to the way variants are reported in the literature

	acceptor_splice = {}
	donor_splice = {}
	cdna2gdna = {}
	cumulative_length = 0
	exon_boundary = 0
	total_length = len(seq)

	for exon in exons:
		acceptor_splice[exon_boundary+1] = exon.acceptor_splice
		cdna2gdna[exon_boundary+1] = exon.start
		cumulative_length += exon.end-exon.start+1

		exon_boundary  = total_length - cumulative_length if reverse else cumulative_length

		donor_splice[exon_boundary] = exon.donor_splice
		cdna2gdna[exon_boundary] = exon.end

	print(cumulative_length, cumulative_length%3, cumulative_length/3)
	biopython_dna = Seq(seq, generic_dna)
	if reverse: biopython_dna = biopython_dna.reverse_complement()
	cdna = str(biopython_dna)
	protein = str(biopython_dna.translate())
	# print(biopython_dna)
	# print("============")
	# print(protein)
	# codon = [cdna[i:i+3] for i in range(0,len(cdna),3)]
	# for i in range(len(protein)):
	# 	print(i, protein[i], codon[i])
	# exit()
	return seq_region_name, cdna, protein, acceptor_splice, donor_splice, cdna2gdna


output_format = '''

abca4_chromosome = "{}"

abca4_cdna = \'\'\'
{}
\'\'\'

abca4_protein =  \'\'\'
{}
\'\'\'

abca4_donor_splice = {}

abca4_acceptor_splice = {}

abca4_cdna2gdna = {}

def get_cdna():
	return abca4_cdna.replace("\\n", "")
	
def get_protein():
	return abca4_protein.replace("\\n", "")
	
def get_codons():
	return [abca4_cdna[i:i+3] for i in range(0,len(abca4_cdna),3)]
	
def get_backward_numbered_acceptor_splice(bdry_pos):
	cumulative_length = max(abca4_acceptor_splice.keys())
	return abca4_acceptor_splice.get(cumulative_length-bdry_pos, None)
	
def get_backward_numbered_donor_splice(bdry_pos):
	cumulative_length = max(abca4_donor_splice.keys())
	return abca4_donor_splice.get(cumulative_length-bdry_pos, None)

'''

########################################
def main():

	print("careful, this will rewrite  utils/abca4_gene.py")
	exit()

	gene = "ABCA4"

	blastdbcmd = "/usr/bin/blastdbcmd"
	cdna_fasta = "/storage/databases/ensembl-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.fa"
	utils_dir = "utils"
	for dep in [blastdbcmd, cdna_fasta, utils_dir]:
		if not os.path.exists(dep):
			print(f"{dep} not found")
			exit()

	db,cursor = abca4_connect()

	# ensmebl id from gene name
	qry = f"select ensembl_gene_id from identifier_maps.hgnc where approved_symbol='{gene}'"
	ensembl_gene_id = hard_landing_search(cursor, qry)[0][0]
	# print(ensembl_gene_id)

	switch_to_db(cursor, "homo_sapiens_core_97_38")

	# coordinates of canonical exons
	# we will just extract positions of exon boundaries on cdna
	[seq_region_id, strand, exons] = get_canonical_exon_bdries(cursor, ensembl_gene_id)

	# canonical cdna and its translation
	[seq_region_name, cdna, protein, acceptor_splice, donor_splice, cdna2gdna] =\
		get_cdna(cursor, blastdbcmd, cdna_fasta, seq_region_id, strand, exons)
	pp = pprint.PrettyPrinter(indent=4)
	with open(f"{utils_dir}/{gene.lower()}_gene.py", "w") as outf:
		newlined_cdna = "\n".join([cdna[i:i+100] for i in range(0,len(cdna),100)])
		newlined_protein = "\n".join([protein[i:i+100] for i in range(0,len(protein),100)])
		outf.write(output_format.format(seq_region_name, newlined_cdna, newlined_protein, pp.pformat(donor_splice),
										pp.pformat(acceptor_splice),  pp.pformat(cdna2gdna)))

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()

