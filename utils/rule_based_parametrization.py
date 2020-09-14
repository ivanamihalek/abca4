
from utils.utils import  *


# aa prp table from http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html
# molmass in Da, columee in cubi \AA,
# hydropathy index from Kyte, J. and Doolittle, R.F., J. Mol. Biol., 157:105-132 (1982).
molecular_mass = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147, 'G': 75, 'H': 155, 'I': 131,
                  'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117}

number_of_atoms = {'A': 13, 'R': 26, 'N': 17, 'D': 16, 'C': 14, 'Q': 20, 'E': 19, 'G': 10, 'H': 20, 'I': 22, 'L': 22,
                   'K': 24, 'M': 20, 'F': 23, 'P': 17, 'S': 14, 'T': 17, 'W': 27, 'Y': 24, 'V': 19}

volume = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5, 'Q': 143.8, 'E': 138.4, 'G': 60.1, 'H': 153.2,
          'I': 166.7, 'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7, 'S': 89.0, 'T': 116.1, 'W': 227.8,
          'Y': 193.6, 'V': 140.0}

hydropathy = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2,
              'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
              'Y': -1.3, 'V': 4.2}

def expression_from_prop_change(aa_from, aa_to):
	[expression, notes]  = [1.0, ""]
	volfrom = volume[aa_from]
	volto = volume[aa_to]
	biggervol = max(volfrom, volto)
	fract = abs(volfrom-volto)/biggervol
	hydro_change = abs(hydropathy[aa_from]-hydropathy[aa_to])
	if hydro_change>3.0:
		expression = 0.25
		notes = "strong expression: large hydropathy change;"
	else:
		if fract<0.2:
			expression = 0.75
			notes = "mild expression: <20% vol change; "
		elif fract<0.3:
			expression = 0.5
			notes = "strong expression:  20%> vol change <30%; "
		else:
			expression = 0.25
			notes = "severe expression:  vol change>30%; "
	return [expression, notes]

def annotate_variable_position(cons_data, protein, protein_domain):
	[expression, transport, notes] = [1.0, 1.0, ""]
	[type_in_specs, substitutions, cons] = cons_data
	[aa_from_long, pos, aa_to_long] = parse_protein(protein)
	aa_from = single_letter_code[aa_from_long]
	aa_to = single_letter_code[aa_to_long]
	if aa_from!=type_in_specs[pos]:
		print(f"type mismatch for pos {pos}: {aa_from} vs {type_in_specs[pos]} in specs")
		exit()
	if aa_to in substitutions[pos]:
		notes = "none: substitution seen"
	else:
		# this isn't really making things any better
		# [expression, notes] = expression_from_prop_change(aa_from, aa_to)
		[expression, notes] = [1.0, ""]

		if protein_domain in ["NBD1", "NBD2"]:
			transport = 0.75
			notes += "mild transport: NBD, novel subs"
		elif protein_domain in ["TMD1", "TMD2"]:
			transport = 0.50
			notes += "strong transport: TMD, novel subs"
		elif protein_domain in ["ECD1", "ECD2"]:
			transport = 0.25
			notes += "strong transport: ECD, novel subs"
		elif protein_domain in ["R"]:
			transport = 0.50
			notes += "severe transport: R, novel subs"
		else:
			transport = 0.75
			notes += f"mild transport: {protein_domain}, novel subs"

	return [expression, transport, notes]


def rule_based_params(variant_info, cons_data):
	[expression, transport, notes] = [1.0, 1.0, ""]

	protein = variant_info["protein"]
	cdna = variant_info["cdna"]
	pd = variant_info["protein_domain"]
	homozygotes =  variant_info["gnomad_homozygotes"]
	conserved_in_para_verts  =  variant_info["conserved_in_para_verts"]
	conserved_in_ortho_verts =  variant_info["conserved_in_ortho_verts"]

	############################################
	# expression
	if "fs" in protein:
		if may_escape_NMD(cdna, protein):
			expression = 0.5
			notes = "strong expression: escapes NMD"
		else:
			expression = 0.0
			notes = "null: fs"
	elif is_splice(cdna, protein):
		expression = 0.0
		notes = "null: near splice"
	elif is_ter(cdna, protein):
		if may_escape_NMD(cdna, protein):
			expression = 0.5
			notes = "strong expression: escapes NMD"
		else:
			expression = 0.0
			notes = "null: early stop"
	elif is_distant_splice(cdna, protein):
		expression = 0.5
		notes = "strong expression: distant splice"
	############################################
	# transport
	elif is_salt_bridge(protein):
		if pd in ["NBD1", "NBD2"]:
			transport = 0.50
			notes = "strong transport: NBD, saltbridge"
		elif pd in ["TMD1", "TMD2"]:
			transport = 0.75 # we don's seem to have a salt bridge here
			notes = "mild transport: TMD, saltbridge"
		elif pd in ["ECD1", "ECD2"]:
			transport = 0.25
			notes = "severe transport: ECD, saltbridge"
		else:
			transport = 0.75 # here neither
			notes = f"mild transport: {pd}, cons in ortho"


	elif is_synonymous(protein):
		notes = "suspicious: synonymous"

	elif in_nucleotide_neighborhood(protein):
		transport = 0.5
		notes = "strong transport: in the nucleotide neighborhood"

	elif conserved_in_para_verts==1:
		if homozygotes and homozygotes>1:
			expression = 0.75
			notes = "mild expression: cons in para but homozygotes exist"
		else:
			expression = 0.5
			notes = "strong expression: cons in para"

	elif conserved_in_ortho_verts==1:
		if pd in ["NBD1", "NBD2"]:
			transport = 0.50
			notes = "strong transport: NBD, cons in ortho"
		elif pd in ["TMD1", "TMD2"]:
			transport = 0.25
			notes = "severe transport: TMD cons in ortho"
		elif pd in ["ECD1", "ECD2"]:
			transport = 0.0
			notes = "no transport: ECD, cons in ortho"
		elif pd in ["R"]:
			transport = 0.25
			notes = "severe transport:  R cons in ortho"
		else:
			transport = 0.50
			notes = f"strong expr and transp: {pd}, cons in ortho"

	elif "ins" in protein or "del" in protein or "dup" in protein:
		if pd=="linker":
			expression = 0.5
			transport = 0.5
			notes = "strong expr and transp: linker region mod"
		else:
			expression = 0.0
			notes = "null: gross structural mod"

	elif "deep" in protein :
		notes = "suspicious: putative deep intronic"

	elif is_missense(protein):
		# for positions that are not completely conserved (not in orthos and not in paras, we took care of that above)
		# check is the subsitution is seen in the alignment or (in which case it might be tolerated)
		# or is  novel - that is, no other sequence has it at this position
		[expression, transport, notes] = annotate_variable_position(cons_data, protein, pd)
	else:
		print(f"warning: failed to annotate {cdna}, {protein} ")

	return [expression, transport, notes]

