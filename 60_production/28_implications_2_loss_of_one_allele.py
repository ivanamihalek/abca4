#!/usr/bin/python3

from utils.simulation import *
import matplotlib.pyplot as plt



#########################################
def main():

	label_format = "%s\na1:  %.0f%%  %.0f%%\na2:  %.0f%%  %.0f%%"
	max_age = 50
	# healthy
	title = "healthy individual"
	params1 = [1.0, 1.0]
	params2 = [1.0, 1.0]
	rpe_baseline = 0.1
	x1, y1 = sim(params1, params2, rpe_baseline, max_age=max_age)
	# shift so it is visible on the graph
	y1_shifted = [y+0.01 for y in y1["rpe"]]
	#variants1 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)

	# no protein product
	title = "no protein product"
	params1 = [0.0, 1.0]
	params2 = [0.0, 1.0]
	rpe_baseline = 0.1
	x2, y2 = sim(params1, params2, rpe_baseline, max_age=max_age)
	#variants2 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)

	title = "loss of product\nfrom one allele"
	params1 = [1.0, 1.0]
	params2 = [0.0, 1.0]
	rpe_baseline = 0.1
	x3, y3 = sim(params1, params2, rpe_baseline, max_age=max_age)
	variants3 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)

	title = "loss of transport\nin one allele"
	params1 = [1.0, 1.0]
	params2 = [1.0, 0.0]
	rpe_baseline = 0.1
	x4, y4 = sim(params1, params2, rpe_baseline, max_age=max_age)
	variants4 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)



	fig, axs = plt.subplots()
	axs.set_ylim(-0.1,1.1)
	axs.text(0.65, 0.68, variants3, transform = axs.transAxes)
	axs.text(0.30, 0.25, variants4, transform = axs.transAxes)
	axs.plot(x1, y1_shifted)
	axs.plot(x2, y2["rpe"])
	axs.plot(x3, y3["rpe"])
	axs.plot(x4, y4["rpe"])



	plt.show()
	return



#########################################
if __name__ == '__main__':
	main()

