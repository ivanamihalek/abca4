
from utils.mysql import *
from math import log10


def human_readble(frequency):
	if not frequency: return "0"
	logten = -log10(frequency)
	if logten>6:
		return "<1:1M"
	if logten>5:
		return "<1:100K"
	if logten>4:
		return "<1:10K"
	if logten>3:
		return "<1:1K"
	if logten>1:
		return "%.0f:1K"%(frequency*1000)
	return "%.0f:10"%(frequency*10)


def variants_from_allele(cursor, allele_id):
	qry = f"select variant_ids from alleles where id = {allele_id}"
	var_ids = [v[1:-1].split("-") for v in hard_landing_search(cursor, qry)[0]]
	for vi in var_ids:
		qry  = "select cdna, protein, gnomad_freq, gnomad_homozygotes, protein_domain, conserved_in_ortho "
		qry += f"from variants where id in ({','.join(vi)})"
		cdna = []
		protein = []
		freqs = []
		homozygs = []
		regions = []
		conserved = []
		for [c, p, f, h, r, cons] in hard_landing_search(cursor, qry):
			cdna.append(f"c.{c}" if c else "np")
			protein.append(f"p.{p}")
			freqs.append(human_readble(f))
			homozygs.append(str(h) if h else "0")
			regions.append(r if r and not "ter" in p.lower() else "none")
			conserved.append("Y" if cons and cons==1 else "n")

	ret = ["; ".join(cdna), "; ".join(protein), "; ".join(freqs)]
	ret.extend(["; ".join(homozygs), "; ".join(regions), ";".join(conserved)])
	return  ret
