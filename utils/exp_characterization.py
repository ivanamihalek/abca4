
characterized_cdna = ["5461-10T>C"]

def characterized_cdna_values(cdna):
	[expression, transport, notes] = [1.0, 1.0, ""]
	if "5461-10T>C" in cdna: # https://pubmed.ncbi.nlm.nih.gov/27775217/
		expression = 0.15
		notes = "severe expression: from exp"

	return [expression, transport, notes]
