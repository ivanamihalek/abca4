
from PIL import Image
import os, re
from utils.annotation import  single_letter_code
from utils.structure import putative_salt_bridge_members, nucleotide_6A_neighborhood
from utils.abca4_gene import abca4_donor_splice

#########################################

def unpack_progression(progression):
	age = []
	va = []
	for agepoint in progression.split(";"):
		[a, v]  = [float(number) for number in agepoint.split(":")]
		age.append(a)
		va.append(v)
	return age, va

#########################################

def parse_protein(protein):
	if not protein: return [None, None, None]
	pattern = re.match('(\D{3})(\d+)(\D{3})', protein.strip().replace("p.",""))
	if not pattern: return [None, None, None]
	aa_from = pattern.group(1).upper()
	pos = int(pattern.group(2))
	aa_to = pattern.group(3).upper()
	return [aa_from, pos, aa_to]

def is_ter(cdna, prot):
	if prot and "ter"in prot.lower(): return True
	return False

def may_escape_NMD(cdna, protein):
	# this just checks 50-55 rule (does not check that the stop codon is actually generated)
	last_splice_position = sorted(abca4_donor_splice.keys())[-1]
	if cdna:
		pattern = re.match('(\d+)\D', cdna.strip().replace("c.",""))
		if not pattern: return False
		pos = int(pattern.group(1))
		if pos>last_splice_position-60: return True
	if protein:
		[aa_from, pos, aa_to]  = parse_protein(protein)
		if pos*3 >last_splice_position-60: return True
	return False

def is_del(cdna, prot):
	if prot and "del"in prot.lower(): return True
	return False

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

def is_salt_bridge(protein):
	if not protein: return False
	pattern = re.match('(\D{3})(\d+)(\D{3})', protein.strip().replace("p.",""))
	if not pattern: return False

	aa_from = pattern.group(1).upper()
	if aa_from.upper()!="CYS": return False

	pos = int(pattern.group(2))
	if not pos in putative_salt_bridge_members: return False

	aa_to = pattern.group(3).upper()
	if aa_to.upper()=="CYS": return False

	return True

def in_nucleotide_neighborhood(protein):
	if not protein: return False
	pattern = re.match('(\D{3})(\d+)(\D{3})', protein.strip().replace("p.",""))
	if not pattern: return False

	pos = int(pattern.group(2))
	if not pos in nucleotide_6A_neighborhood: return False

	aa_from = pattern.group(1).upper()
	aa_to = pattern.group(3).upper()
	if aa_from==aa_to:  return False

	return True



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


def panic(panic_args):
	print(panic_args)
	exit()

