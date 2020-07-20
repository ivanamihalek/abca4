
from PIL import Image
import os, re
from utils.annotation import  single_letter_code

#########################################

def is_ter(cdna, prot):
	if prot and "ter"in prot.lower(): return True
	return False

def is_del(cdna, prot):
	if prot and "del"in prot.lower(): return True
	return False

def parse_protein(protein):
	if not protein: return [None, None, None]
	pattern = re.match('(\D{3})(\d+)(\D{3})', protein.strip().replace("p.",""))
	if not pattern: return [None, None, None]
	aa_from = pattern.group(1).upper()
	pos = int(pattern.group(2))
	aa_to = pattern.group(3).upper()
	return [aa_from, pos, aa_to]

def is_missense(protein):
	[aa_from, pos, aa_to] = parse_protein(protein)
	if not aa_from or not aa_to: return False
	if aa_from!=aa_to: return True
	return False

def is_synonymous(protein):
	[aa_from, pos, aa_to] = parse_protein(protein)
	if not aa_from or not aa_to: return False
	if aa_from==aa_to: return True
	return False

# this is not very reliable:
# we have a case wiht two "misfolder" alleles that is over 40 at onset
def is_misfolder(cdna, protein_vars):
	if not protein_vars: return False
	if "_" in protein_vars: return False
	for protein in protein_vars.split(";"):
		pattern = re.findall('(\D{3})(\d+)(\D{3})', protein.strip().replace("p.",""))
		if not pattern: return False
		for match in pattern:
			aa_from = match[0].upper()
			pos = match[1]
			aa_to = match[2].upper()
			# print(f"{protein}  {aa_from}  {aa_to} ")
			# print(f"{single_letter_code[aa_from]}{pos}{single_letter_code[aa_to]}")
			if f"{single_letter_code[aa_from]}{pos}{single_letter_code[aa_to]}" in ["A1357T", "A1794P", "L2027F", "R2077W"]:
				return True
	return False


def is_exotic(cdna, protein):
	if not cdna: return False
	if "5461-10" in cdna: return True
	return False


def parse_splice(cdna):
	pattern = re.findall('(\d+)[\-+](\d+)(\D)', cdna)
	min_match = 100
	if pattern:
		min_match = 50
		for match in pattern:
			if int(match[1])<min_match: min_match=int(match[1])
	return min_match

def is_splice(cdna, protein):
	if not cdna: return False
	splice_position = parse_splice(cdna)
	if splice_position<3: return True

	# last nucleotide before splice
	pattern = re.findall('(\D+)(\d+)(\D+)splice', protein)
	if pattern: return True

	return False

def is_distant_splice(cdna, protein):
	if not cdna: return False
	splice_position = parse_splice(cdna)
	if splice_position<13: return True

	return False


# if "del"in prot.lower(): return True
# if "5461-10T>C" in cdna: return True # we have exp proof of this one
def is_null(cdna, prot, filter):
	if len(filter)==0: return True
	if "ter" in filter and is_ter(cdna, prot): return True
	if "del" in filter and is_del(cdna, prot): return True
	if "splice" in filter and is_splice(cdna, prot): return True
	if "exotic" in filter and is_exotic(cdna, prot): return True
	if "misfolder" in filter and is_misfolder(cdna, prot): return True
	return False


def convert_to_jpg(fnm):
	# note e includes "."
	f, e = os.path.splitext(fnm)
	new_fnm = f + ".jpg"
	try:
		Image.open(fnm).convert('RGB').save(new_fnm)
	except IOError as err:
		print("cannot convert", fnm, err)
		return None
	os.remove(fnm)
	return


def convert_frames(frames_home, subdir):
	outdir = "{}/{}".format(frames_home, subdir)
	for path, subdirs, files in os.walk(outdir):
		print(files[:5])
		for fnm in files:
			if fnm[-3:]=="png":
				convert_to_jpg("{}/{}".format(path,fnm))
	return


def subdir_prep(frames_home, subdir):
	outdir = "{}/{}".format(frames_home, subdir)
	#if os.path.exists(outdir): shutil.rmtree(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

