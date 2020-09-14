#!/usr/bin/python3

from utils.simulation import *
import matplotlib.pyplot as plt



#########################################
def main():

	label_format = "%s\na1:  %.0f%%  %.0f%%\na2:  %.0f%%  %.0f%%"
	label_format_intervention = "%s\na1:  %.0f%%  %.0f%%\na2:  %.0f%%  %.0f%%\nintervention:  %.0f%%  %.0f%%"
	max_age = 50
	# healthy
	title = "healthy individual"
	params1 = [1.0, 1.0]
	params2 = [1.0, 1.0]
	rpe_baseline = 0.1
	x1, y1 = sim(params1, params2, rpe_baseline, max_age=max_age)
	# shift so it is visible on the graph
	y1_shifted = [y+0.01 for y in y1["rpe"]]
	variants1 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)

	# no protein product
	title = "no protein product"
	params1 = [0.0, 1.0]
	params2 = [0.0, 1.0]
	rpe_baseline = 0.1
	x2, y2 = sim(params1, params2, rpe_baseline, max_age=max_age)
	#variants2 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)


	title = "loss of transport\nin one allele"
	params1 = [1.0, 1.0]
	params2 = [1.0, 0.0]
	rpe_baseline = 0.1
	x3, y3 = sim(params1, params2, rpe_baseline, max_age=max_age)
	variants3 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)
	color3 =  '#c99797'  # super pale red

	# intervention
	intervention_age = 5
	params1 = [1.0, 1.0]
	params2 = [1.0, 0.0]
	params3 = [1.0, 1.0]
	rpe_baseline = 0.1
	x4, y4 = sim_intervention(params1, params2,  params3, rpe_baseline, max_age=max_age, intervention_age=intervention_age)
	variants4 = f"intervention at age {intervention_age}"
	color4 = '#e67778' # pale red

	# intervention
	intervention_age = 2
	params1 = [1.0, 1.0]
	params2 = [1.0, 0.0]
	params3 = [1.0, 1.0]
	rpe_baseline = 0.1
	x5, y5 = sim_intervention(params1, params2,  params3, rpe_baseline, max_age=max_age, intervention_age=intervention_age)
	color5 = '#d62728'  # default red

	# intervention
	intervention_age = 5
	params1 = [1.0, 1.0]
	params2 = [1.0, 0.0]
	params3 = [1.0, 1.0]
	silenced_alleles = [1]
	rpe_baseline = 0.1
	x6, y6 = sim_intervention(params1, params2,  params3, rpe_baseline,
	                          max_age=max_age, intervention_age=intervention_age, silenced_alleles=silenced_alleles)
	color6 = '#e67778' # pale red


	# intervention
	intervention_age = 2
	params1 = [1.0, 1.0]
	params2 = [1.0, 0.0]
	params3 = [1.0, 1.0]
	silenced_alleles = [1]
	rpe_baseline = 0.1
	x7, y7 = sim_intervention(params1, params2,  params3, rpe_baseline,
	                          max_age=max_age, intervention_age=intervention_age, silenced_alleles=silenced_alleles)
	color7 = '#d62728'  # default red


	fig, axs = plt.subplots()
	axs.set_ylim(-0.1,1.1)
	axs.plot(x1, y1_shifted)
	axs.plot(x2, y2["rpe"])
	axs.text(0.30, 0.25, variants4, transform = axs.transAxes)
	axs.plot(x3, y3["rpe"], color=color3, linestyle='dashed')
	axs.plot(x4, y4["rpe"], color=color4, linestyle='dashed')
	axs.plot(x5, y5["rpe"], color=color5, linestyle='dashed')
	axs.plot(x6, y6["rpe"], color=color6)
	axs.plot(x7, y7["rpe"], color=color7)


	plt.show()
	return



#########################################
if __name__ == '__main__':
	main()

