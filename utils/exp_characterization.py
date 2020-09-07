
characterized_protein = ['Leu2027Phe', 'Arg2030Gln', 'Arg602Gln', 'Asp1532Tyr']

def characterized_prot_values(protein):
	[expression, transport, notes] = [1.0, 1.0, ""]
	# todo move to database
	# todo make systeatic: how do measured values contribute to this estimate

	if 'Gly72Arg' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		transport = 0.1
		notes = "severe transport: from exp"

	elif 'Met448Lys' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		transport = 0.4
		notes = "strong transport: from exp"

	elif 'Leu541Pro' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		transport = 0.1
		notes = "severe transport: from exp"

	elif 'Val552Ile' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		transport = 0.5
		notes = "strong transport: from exp"

	elif 'Asn965Ser' in protein: # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5886264/
		expression = 0.5
		transport = 0.4
		notes = "strong transport: from exp"

	elif 'Ala1038Val' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		expression = 0.75
		transport = 0.4
		notes = "strong transport: from exp"

	elif 'Gly1091Asp' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		expression = 0.75
		transport = 0.75
		notes = "mild transport: from exp"

	elif 'Ala1357Thr' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		expression = 0.5
		transport = 0.75
		notes = "strong expression: from exp"

	elif 'Ala1794Pro' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		expression = 0.5
		transport = 0.4
		notes = "strong transport: from exp"

	elif 'Gly1961Asp' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		transport = 0.4
		notes = "strong transport: from exp"

	elif 'Leu2027Phe' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		expression = 0.5
		transport = 0.4
		notes = "strong transport: from exp"

	elif 'Arg2077Trp' in protein:  # https://iovs.arvojournals.org/article.aspx?articleid=2680974
		expression = 0.5
		transport = 0.4
		notes = "strong transport: from exp"



	# elif 'Leu2027Phe' in protein: # I cannot find the reference for this
	# 	# the experiment says this is a misfolder, but does not seem to be
	# 	expression = 0.75
	# 	notes = "mild expression: misfolder"
	# elif 'Arg2030Gln' in protein:
	# 	expression = 0.75
	# 	notes = "mild expression: from exp"
	# elif 'Arg602Gln' in protein:
	# 	# 602 is and ECD1 position conserved in vertebrates
	# 	# accordingly  Arg602Trp is assgned no-transport value
	# 	# which corresponds t the observed onset range og 6-26 yrs in carrrier
	# 	# however Arg602Gln carriers  have onsets of 31, 38, 57
	# 	# the 31yr onset is a homozygote in this variant
	# 	transport = 0.50
	# 	notes = "strong transport: ECD1, from exp"
	# elif 'Asp1532Tyr' in protein:
	# 	# 1532 is and ECD2 position conserved in vertebrates
	# 	# however Asp1532Tyr has only two cases, with onset  at 36 and 37
	# 	transport = 0.50
	# 	notes = "strong transport: ECD2, from exp"

	return [expression, transport, notes]

characterized_cdna = ["5461-10T>C"]

def characterized_cdna_values(cdna):
	[expression, transport, notes] = [1.0, 1.0, ""]
	if "5461-10T>C" in cdna: # https://pubmed.ncbi.nlm.nih.gov/27775217/
		expression = 0.15
		notes = "severe expression: from exp"

	return [expression, transport, notes]
