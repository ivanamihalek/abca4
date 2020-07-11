
from utils.mysql import *

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
