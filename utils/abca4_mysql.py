
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
		qry = f"select cdna, protein, gnomad_freq, gnomad_homozygotes from variants where id in ({','.join(vi)})"
		cdna = []
		protein = []
		freqs = []
		homozygs = []
		for [c, p, f, h] in hard_landing_search(cursor, qry):
			cdna.append(f"c.{c}" if c else "np")
			protein.append(f"p.{p}")
			freqs.append(human_readble(f))
			homozygs.append(str(h) if h else "0")
	return "; ".join(cdna), "; ".join(protein), "; ".join(freqs), "; ".join(homozygs)
