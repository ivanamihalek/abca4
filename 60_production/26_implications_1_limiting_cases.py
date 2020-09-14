#!/usr/bin/python3

from utils.simulation import *
import matplotlib.pyplot as plt



#########################################
def main():

	label_format = "%s\na1:  %.0f%%  %.0f%%\na2:  %.0f%%  %.0f%%"

	# healthy
	title = "healthy individual"
	params1 = [1.0, 1.0]
	params2 = [1.0, 1.0]
	rpe_baseline = 0.1
	x1, y1 = sim(params1, params2, rpe_baseline, max_age=200)
	variants1 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)


	# no protein product
	title = "no protein product"
	params1 = [0.0, 1.0]
	params2 = [0.0, 1.0]
	rpe_baseline = 0.1
	x2, y2 = sim(params1, params2, rpe_baseline, max_age=200)
	variants2 = label_format % (title, params1[0]*100, params1[1]*100, params2[0]*100, params2[1]*100)


	fig, axs = plt.subplots()
	axs.set_ylim(-0.1,1.1)
	axs.text(0.6, 0.75, variants1, transform = axs.transAxes)
	axs.text(0.15, 0.25, variants2, transform = axs.transAxes)
	axs.plot(x1, y1["rpe"])
	axs.plot(x2, y2["rpe"])



	plt.show()
	return



#########################################
if __name__ == '__main__':
	main()

