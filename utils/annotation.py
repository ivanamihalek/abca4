
import re
from math import floor

from Bio.Seq  import Seq, MutableSeq
# mutable seq is not particulalry useful for deleations because it can only delete a single nucleotide
# biopython_dna = MutableSeq(seq, generic_dna)
# biopython_dna.pop(position)
from Bio.Alphabet import generic_dna, generic_protein
from Bio import BiopythonWarning
import warnings
warnings.simplefilter('ignore', BiopythonWarning)

from utils.abca4_gene import  *

single_letter_code = {'GLY':'G', 'ALA':'A',  'VAL':'V', 'LEU':'L', 'ILE':'I',
		 'MET': 'M', 'PRO': 'P',  'TRP': 'W', 'PHE':'F', 'SER':'S',
		 'CYS': 'C', 'THR': 'T',  'ASN': 'N', 'GLN':'Q', 'TYR':'Y',
		 'LYS': 'K', 'ARG': 'R',  'HIS': 'H', 'ASP':'D', 'GLU':'E','TER':'*'}

three_letter_code = dict([(single_letter_code[k],k.capitalize()) for k in single_letter_code.keys()])
# note that we have 3 different glyphs in use for the minus character
sign_regex = r'[\-+−–]'
minus_character = "-−–"


def indel(seq, cdna_variant, original_protein, verbose=False):
	# return
	# trying to accomodate various non-standard ways of writing down the variant
	if "indel" in cdna_variant:
		field = cdna_variant.split("indel")
	elif "delins" in cdna_variant:
		# this is actually the recommended term http://varnomen.hgvs.org/recommendations/DNA/variant/delins/
		field = cdna_variant.split("delins")
	else:
		print(f"what the heck is {cdna_variant}?")
		exit()
	if "_" in cdna_variant:
		position = [int(i) - 1 for i in field[0].split("_")]
	else:
		position = [int(field[0]) - 1, int(field[0])- 1]

	inserted_nt = field[1]
	if len(inserted_nt)==0:
		print(f"no sequence in indel {cdna_variant}?")
		exit()


	insert_length = len(inserted_nt)
	if insert_length%3==0:
		inserted_aa = str(Seq(inserted_nt, generic_dna))
		protein_effect = f"inserted {inserted_aa}"
	else:
		biopython_dna = Seq(seq[:position[0]]+inserted_nt+seq[position[1]+1:], generic_dna)
		new_protein = str(biopython_dna.translate())
		first = None
		fs_start = 1
		fs_end = len(original_protein)
		for i in range(len(new_protein)):
			if original_protein[i] != new_protein[i]:
				if first is None:
					first= f"{three_letter_code[original_protein[i]]}{i+1}{three_letter_code[new_protein[i]]}"
					fs_start = i
			if new_protein[i]=="*" or original_protein[i]=="*":
				fs_end = i
				break
		protein_effect = f"{first}fsTer{fs_end-fs_start+1}"
	if verbose: print(">>>>>>>>>>>>>>>>", position, inserted_nt, protein_effect)
	return protein_effect


def insert(seq, cdna_variant, original_protein, verbose=False):
	# return
	# trying to accomodate various non-standard ways of writing down the variant
	field = cdna_variant.split("ins")
	if "_" in cdna_variant:
		# example 3211_3212insGT
		position = [int(i) - 1 for i in field[0].split("_")]
	else:
		# c.3211insGT
		position  = [int(field[0]) - 1, int(field[0]) - 1]

	inserted_nt = field[1]

	insert_length = len(inserted_nt)
	if insert_length%3==0:
		inserted_aa = str(Seq(inserted_nt, generic_dna).translate())
		protein_effect = f"inserted {inserted_aa}"
	else:
		biopython_dna = Seq(seq[:position[0]+1]+inserted_nt+seq[position[1]:], generic_dna)
		new_protein = str(biopython_dna.translate())
		first = None
		fs_start = 1
		fs_end = len(original_protein)
		for i in range(len(new_protein)):
			if original_protein[i] != new_protein[i]:
				if first is None:
					first= f"{three_letter_code[original_protein[i]]}{i+1}{three_letter_code[new_protein[i]]}"
					fs_start = i
			if new_protein[i]=="*" or original_protein[i]=="*":
				fs_end = i
				break
		protein_effect = f"{first}fsTer{fs_end-fs_start+1}"
	if verbose: print(">>>>>>>>>>>>>>>>", protein_effect)
	return protein_effect


def deletion(seq, cdna_variant, original_protein, verbose=False):	#return
	# trying to accomodate various non-standard ways of writing down the variant
	if "_" in cdna_variant:
		# example 4663_4664del
		# however, can also have 850_857delATTCAAGA
		if verbose: print(cdna_variant)
		position = [int(i) - 1 for i in cdna_variant.split("del")[0].split("_")]
	else:
		# I've seen 4739del, but also 4739delT for the same thing
		position  = [int(cdna_variant.split("del")[0]) - 1]

	if len(position)==1: position.append(position[0])
	deletion_length = position[1] - position[0] + 1
	# did we delete a pece of protein seq without frameshift?
	# we are in offset 0  system - if the deletion starts at position[0]==3, means there are 0, 1, and 2 left behind
	if position[0]%3==0 and (len(seq)-position[1]-1)%3==0:
		del_start_aa_type = three_letter_code[str(Seq(seq[position[0]:position[0]+3], generic_dna).translate())]
		del_end_aa_type   = three_letter_code[str(Seq(seq[position[1]-2:position[1]+1], generic_dna).translate())]
		protein_effect = f"{del_start_aa_type}{int(position[0]/3)+1}_{del_end_aa_type}{int(position[1]/3)+1}del"
	else:
		biopython_dna = Seq(seq[:position[0]]+seq[position[1]+1:], generic_dna)
		new_protein = str(biopython_dna.translate())
		first = None
		fs_start = 1
		fs_end = len(original_protein)
		for i in range(len(new_protein)):
			if original_protein[i] != new_protein[i]:
				if first is None:
					first= f"{three_letter_code[original_protein[i]]}{i+1}{three_letter_code[new_protein[i]]}"
					fs_start = i
			if new_protein[i]=="*" or original_protein[i]=="*":
				fs_end = i
				break
		protein_effect = f"{first}fsTer{fs_end-fs_start+1}"
	if verbose: print(">>>>>>>>>>>>>>>>", protein_effect)
	return protein_effect


def duplication(seq, cdna_variant, original_protein, verbose=False):
	if "_" in cdna_variant:
		# example 2169_2172dup
		# split in case somebody gave us the actual duplicated sequence
		if verbose: print(cdna_variant)
		position = [int(i)  for i in cdna_variant.split("dup")[0].split("_")]
		cdna_variant_as_insert = f"{position[1]}ins{seq[position[0]-1:position[1]]}"
	else:
		# 4604dup
		position  = int(cdna_variant.split("dup")[0])
		cdna_variant_as_insert = f"{position}ins{seq[position-1]}"

	return insert(seq, cdna_variant_as_insert, original_protein, verbose)



def point_mutation(seq, cdna_variant, original_protein, verbose=False):
	#return ""
	if not cdna_variant: return ""
	variant_parse = re.match(r'(\d+)([ACTG])>([ACTG])', cdna_variant)
	if not variant_parse:
		print("not parseable:", cdna_variant)
		exit()
	pos = int(variant_parse.group(1))
	nt_from = variant_parse.group(2)
	nt_to = variant_parse.group(3)
	pos -= 1  # 0-offset
	# sanity
	if seq[pos]!=nt_from:
		if verbose: print(f"seq mismatch in from: {nt_from} should be {seq[pos]}")
	protein_pos = floor(pos/3)
	biopython_dna = MutableSeq(seq)
	biopython_dna[pos] = nt_to
	new_protein = str(biopython_dna.toseq().translate())

	protein_effect = f"{three_letter_code[original_protein[protein_pos]]}{protein_pos+1}{three_letter_code[new_protein[protein_pos]]}"
	# special warning: if the last G before splice is mutated, that might affect splicing
	if pos+1 in abca4_donor_splice and nt_from=="G": protein_effect += "splice"


	return protein_effect

###################################
def splice_site_IVS(cdna_variant, verbose=False):
	# check only that the exon boundary is recognizable
	if verbose: print("\n ====> ", cdna_variant)
	# parse the variant
	variant_parse = re.match(r'IVS(\d+)([\-+−–])(\d+)([ACTG])>([ACTG])', cdna_variant)
	if not variant_parse:
		print("not parseable IVS:", cdna_variant)
		exit()
	exon_number = int(variant_parse.group(1))
	direction = variant_parse.group(2)
	pos       = int(variant_parse.group(3))
	nt_from   = variant_parse.group(4)
	nt_to     = variant_parse.group(5)

	ss_length = len(list(abca4_acceptor_splice.values())[0]) # they should all be the same length
	if direction in minus_character:
		splice = abca4_acceptor_splice
		pos = ss_length-pos
	else:
		splice = abca4_donor_splice
		pos -= 1  # 0-offset
		exon_number -= 1

	cdna_position_of_the_intron = sorted(splice.keys())[exon_number]
	splice_seq = splice[cdna_position_of_the_intron]
	if splice_seq[pos] == nt_from:
		cdna_variant = f"{cdna_position_of_the_intron}{direction}{pos+1}{nt_from}{nt_to}"
	else:
		cdna_variant = "splice   not reproducible"
		if verbose:
			print(exon_number, direction, pos, nt_from, nt_to)
			print(cdna_position_of_the_intron)
			print(splice_seq[:pos], splice_seq[pos], splice_seq[pos+1:])
			print(Seq(splice_seq).reverse_complement())

	return cdna_variant

###################################
def splice_site_del(cdna_variant, verbose=False):
	# check only that the exon boundary is recognizable
	if verbose: print("\n ====> ", cdna_variant)
	# parse the variant
	variant_parse = re.match(r'(\d+)([\-+−–])(\d+)_(\d+)([\-+−–])(\d+)del', cdna_variant)
	if variant_parse:
		exon_bdry_1  = int(variant_parse.group(1))
		direction_1  = variant_parse.group(2)
		pos_1        = int(variant_parse.group(3))
		exon_bdry_2  = int(variant_parse.group(4))
		direction_2  = variant_parse.group(5)
		pos_2        = int(variant_parse.group(6))
		if exon_bdry_1!=exon_bdry_2 or direction_1!=direction_2:
			print("unexepected format for deletion in splice site:", cdna_variant)
			exit()
	else:
		variant_parse = re.match(r'(\d+)([\-+−–])(\d+)del', cdna_variant)
		if variant_parse:
			exon_bdry_1  = int(variant_parse.group(1))
			direction_1  = variant_parse.group(2)
			pos_1        = int(variant_parse.group(3))
		else:
			print("not parseable deletion in splice site:", cdna_variant)
			exit()

	ss_length = len(list(abca4_acceptor_splice.values())[0]) # they should all be the same length
	if direction_1 in minus_character:
		splice = abca4_acceptor_splice
		pos_1 = ss_length-pos_1
	else:
		splice = abca4_donor_splice

	effect = "splice"
	# if the variant is further away than the sequence we have stored return "deep intronic"
	if pos_1<0 or pos_1>=ss_length:
		effect = "splice   corrected: deep intronic"
		if verbose: print("      ================>", effect)
		splice_seq = splice.get(exon_bdry_1, "not found")
		if splice_seq=="not found":
			effect = "splice   not reproducible"

	return effect


###################################
def splice_site(cdna_variant, verbose=False):
	if verbose: print("\n ====> ", cdna_variant)
	# parse the variant
	variant_parse = re.match(r'(\d+)([\-+−–])(\d+)([ACTG])>([ACTG])', cdna_variant)
	if not variant_parse:
		print("not parseable:", cdna_variant)
		exit()
	exon_bdry = int(variant_parse.group(1))
	direction = variant_parse.group(2)
	pos       = int(variant_parse.group(3))
	nt_from   = variant_parse.group(4)
	nt_to     = variant_parse.group(5)

	if verbose: print(" ====> ", cdna_variant, exon_bdry, direction, pos, nt_from, nt_to)
	ss_length = len(list(abca4_acceptor_splice.values())[0]) # they should all be the same length
	if direction in minus_character:
		splice = abca4_acceptor_splice
		pos = ss_length-pos
	else:
		pos -= 1  # 0-offset
		splice = abca4_donor_splice

	effect = "splice"
	# if the variant is further away than the sequence we have stored return "deep intronic"
	if pos<0 or pos>=ss_length:
		effect = "splice   corrected: deep intronic"
		if verbose: print("      ================>", effect)
	else:
		splice_seq = splice.get(exon_bdry, "not found")
		if splice_seq=="not found":
			# if exon_bdry not found see if we are counting the exons backwards
			if verbose: print("      ================> exon bdry not found - try reverse numbering of exons")
			if direction in minus_character:
				splice_seq = get_backward_numbered_acceptor_splice(exon_bdry)
			else:
				splice_seq = get_backward_numbered_donor_splice(pos)
			if splice_seq[pos] == nt_from:
				if verbose: print("      ====> the reverse numbering of exons seems to work")
				if verbose: print(splice_seq,  splice_seq[pos])
				# TODO: what shoudl be the actual splice number here
				effect = "splice   corrected: "
			else:
				# if there is nucleotdie mismatch, see if we are counting the exons backwards
				if verbose: print("      ============================>", splice_seq,  splice_seq[pos])
				if verbose: print("      ============================> nt mismatch - try reverse complement of the sequence")
				splice_seq = str(Seq(splice_seq).reverse_complement())
				if splice_seq[pos] == nt_from:
					if verbose: print("      ============> reverse complement seems to work")
					if verbose: print(splice_seq,  splice_seq[pos])
					# TODO: what shoudl be the actual  variant here
					effect = "splice   corrected: "
				else:
					if verbose: print("      ============> reverse complement of the sequence did not work either")
					effect = "not reproducible"
		else:
			if splice_seq[pos] == nt_from:
				if verbose: print(splice_seq,  splice_seq[pos])
			else:
				# if there is nucleotdie mismatch, see if we are counting the exons backwards
				if verbose: print("      ============================>", splice_seq,  splice_seq[pos])
				if verbose: print("      ============================> nt mismatch - try reverse complement of the sequence")
				splice_seq = str(Seq(splice_seq).reverse_complement())
				if splice_seq[pos] == nt_from:
					if verbose: print("      ============> reverse complement seems to work")
					if verbose: print(splice_seq,  splice_seq[pos])
					# TODO: what shoudl be the actual  variant here
					effect = "splice   corrected: "
				else:
					if verbose: print("      ============> reverse complement of the sequence did not work either")
					effect = "splice   not reproducible"

	return effect


def mutation_effect(cdna_variant):

	cdna_variant = cdna_variant.replace(" ", "")
	# note we have 3 different glyphs in use for the minus character (see def for sign_regex on the top of the file)
	if re.findall(sign_regex, cdna_variant):
		if "del" in cdna_variant:
			return splice_site_del(cdna_variant)
		if "IVS" in cdna_variant:
			return splice_site_IVS(cdna_variant)
		else:
			return splice_site(cdna_variant)

	# everything else except splice site
	seq = get_cdna()
	original_protein = str(Seq(seq).translate())

	if "indel" in cdna_variant or "delins" in cdna_variant:
		return indel(seq, cdna_variant, original_protein)

	if "del" in cdna_variant:
		return deletion(seq, cdna_variant, original_protein)

	if "ins" in cdna_variant:
		return insert(seq, cdna_variant, original_protein)

	if "dup" in cdna_variant:
		return duplication(seq, cdna_variant, original_protein)

	if ">" in cdna_variant:
		return point_mutation(seq, cdna_variant, original_protein)

	return ""

