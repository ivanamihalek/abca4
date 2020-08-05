#! /usr/bin/python3
from math import exp
import numpy as np
import matplotlib.pyplot as plt

'''
"Populations" here refers to the populations of ABCA4 alleles present.
We take there are two, and their ratio depends on the gene's ability to produce a functioning protein.


1) First, the model asumes that the total throughput of ATR through the membrane is the weighted sum of throughput 
   effciencies of each allele. It might be affected by the mutations in the TM chanel, NBD or ECD.
   
   total_throughput = \sum fraction_i x transport_efficiency_i

   The ability of each mutant protein to transport ATR is independent of the rate at which it is delivered
   to the membrane. For example, a splice mutant may have hard time stitching itself together, 
   but once it does, it transports just as well as the wild type.  
 

1) Second, the capability of an ABCA4 allele to express, splice, and fold 
   combines into the rate at which its product is delivered and incorporated into disc membrane of the rod cells. 
   The relative abundance (relative to each other; fraction)i, above) of products of the two alleles incorporated 
   in the disc membrane   is proportional to their relative rate of delivery.  
   In the calling method we see something like this:
	
	delivery_rate = [alpha[0]*exp(-age/beta[0]), alpha[1]*exp(-age/beta[1])]
	
	Alpha combines the ability to splice and fold, whereas the overall power to express, encoded in beta, diminshes
	with the dying RPE layer. Beta might actually be a function of the background of each patient. Which leads us to  ..
	
   
3) Third, the ability of the RPE layer to survive is directly proportional to the combined ATR throughput of 
   all ABCA4 units incorporated in the discs. In a feedback fashion, the rate of delivery of both 
   alleles is diminished in time, proportional to the fraction of RPE dying out.

'''

alpha_wt = 0.05
# beta  =  RPE decay parameter
beta_wt  = 200

max_steps = 1000
destroyable_fraction = 0.5

def progression_sim(max_age, alpha_fraction, transport_efficiency, rpe_baseline = 0.1, verbose=False):

	x = []
	y = {"throughput":[], "fraction_0":[], "fraction_1":[] , "rpe":[]}
	beta = beta_wt
	alpha = [alpha_wt*alpha_fraction[0], alpha_wt*alpha_fraction[1]]

	for age in range(max_age):
		# we shouldn't divide by 0
		f = exp(-(age/beta)**2)
		rpe =2*f/(1+f)
		#rpe = f # doesn't really make much difference - the 2*f/(1+f) is a bit "rounder"

		delivery_rate = [alpha[0]*rpe, alpha[1]*rpe]
		fraction = sim_core(delivery_rate, plot=False, verbose=verbose)
		throughput = fraction[0]*transport_efficiency[0] + fraction[1]*transport_efficiency[1]
		beta *= ((1-destroyable_fraction) + destroyable_fraction*throughput)

		if verbose: print("%2d"%age, ["%.2f"%f for f in fraction], "throughput: %.2f"%throughput)
		x.append(age)
		y["throughput"].append(throughput)
		y["fraction_0"].append(fraction[0])
		y["fraction_1"].append(fraction[1])
		# rescale rpe to have some baseline
		y["rpe"].append(rpe_baseline + (1-rpe_baseline)*rpe)
		#y["rpe"].append(rpe) # this just does not reporoduse the observed behavior

	return x, y


def sim_core (delivery_rate, plot=False, verbose=False):

	number_of_populations = len(delivery_rate)
	# relative fraction of each population in the membrane
	# for exaoke, 0.3 is of type 1, 0.4 is of type 2, and 0.3 is unoccupied
	incorporated_fraction = [0.0]*number_of_populations
	delta = [1.0]*number_of_populations
	points = []
	not_incorporated = [0.0]*number_of_populations
	step = 0
	total_fraction_occupied  = 0

	index_iterator = range(number_of_populations)

	while sum(delta)>1.e-6 and step<max_steps:
		# the probabilty of being incorporated in the membrane is proportional to the size of the available pool
		available_pool = [not_incorporated[i] + delivery_rate[i] for i in index_iterator]
		# and the probability that the space in the mebrane is available
		free_space   = (1.0 - total_fraction_occupied)
		# increase of each  population in the membrane  - this goes to zero as the free space does
		delta = list(map(lambda a: a*free_space, available_pool))
		incorporated_fraction = [incorporated_fraction[i] + delta[i] for i in index_iterator]
		total_fraction_occupied = sum(incorporated_fraction)

		points.append(incorporated_fraction.copy())
		step +=1

	if verbose:
		print(f"loop exited after {max_steps}")
		print("sum delta %.1e  total frac occupied %.3f"%(sum(delta), total_fraction_occupied))
		f1 = points[-1][0]
		f2 = points[-1][1]
		if f1+f2> 0.001:
			print("relative contribution of f1 %.3f"%(f1/(f1+f2)))
			print("relative contribution of f2 %.3f"%(f2/(f1+f2)))
		d1 = delivery_rate[0]
		d2 = delivery_rate[1]
		if d1+d2> 0.001:
			print("d1 relative%.3f"%(d1/(d1+d2)))
			print("d2 relative%.3f"%(d2/(d1+d2)))
		print()
	if plot:
		x = np.linspace(0, 1, len(points))
		y0 = [f[0] for f in points]
		y1 = [f[1] for f in points]
		sumy = [f[0]+f[1] for f in points]
		fig, ax = plt.subplots()
		line1, = ax.plot(x, y0, label='f1')
		line2, = ax.plot(x, y1, label='f2')
		line2, = ax.plot(x, sumy, label='sum')
		ax.legend()
		plt.show()

	return points[-1]

def main():


	delivery_rate = [0.3, 0.5]
	for i in range(10):
		delivery_rate = [d/2 for d in delivery_rate]
		sim_core (delivery_rate, plot=False, verbose=True)

	return

if __name__ == "__main__":
	main()