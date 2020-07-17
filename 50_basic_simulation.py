#! /usr/bin/python3

from utils.simulation import progression_sim
import matplotlib.pyplot as plt


def main ():


	# the values for one case/patient
	alpha_fraction = [0.3, 0.3]
	transport_efficiency = [1.0, 0.5]
	max_age = 80
	x, y  = progression_sim(max_age, alpha_fraction, transport_efficiency, verbose=False)

	fig, ax = plt.subplots()
	plt.ylim(0,1.1)
	ax.set_title("alpha = [{}, {}], e = [{}, {}]".
				format( "%.3f"%alpha_fraction[0], "%.3f"%alpha_fraction[1],
	                    "%.2f"%transport_efficiency[0], "%.2f"%transport_efficiency[1]))
	line1, = ax.plot(x, y["throughput"], label='throughput')
	line2, = ax.plot(x, y["fraction_0"], dashes=[6, 2], label='f1')
	line3, = ax.plot(x, y["fraction_1"], dashes=[6, 2], label='f2')
	line4, = ax.plot(x, y["rpe"],  label='RPE')
	ax.legend()
	plt.show()
	return


####################
if __name__ == "__main__":
	main()
