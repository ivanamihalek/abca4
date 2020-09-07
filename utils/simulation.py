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
	y = {"throughput":[], "fraction_0":[], "fraction_1":[], "rpe":[]}
	beta = beta_wt
	alpha = [alpha_wt*alpha_fraction[0], alpha_wt*alpha_fraction[1]]

	for age in range(max_age):
		# we shouldn't divide by 0
		f = exp(-(age/beta)**2)
		#rpe =2*f/(1+f)
		rpe = f # doesn't really make much difference - the 2*f/(1+f) is a bit "rounder"

		delivery_rate = [alpha[0]*rpe, alpha[1]*rpe]
		fraction = sim_core(delivery_rate, plot=False, verbose=verbose)
		# there is something wrong wtih this  - somteims it misses spectacularly
		# it looks in particular in the cases when a throughput of a variant
		# has to be suppresed all the way to zero (unless there is a bug somewhere)
		# fraction = filling_fractions_for_const_delivery_rates(delivery_rate)
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

########################################
def sim(var1_param, var2_param, rpe_baseline, max_age):
	# the values for one case/patient
	# alpha is the abilty to fold and incorporate
	# "fraction" refers to the fraction of the wild-type capability
	expression_fraction = [var1_param[0], var2_param[0]]
	transport_fraction = [var1_param[1], var2_param[1]]
	x, y  = progression_sim(max_age, expression_fraction, transport_fraction, rpe_baseline=rpe_baseline, verbose=False)
	return x, y

########################################
def filling_fractions_for_const_delivery_rates(delivery_rate):
	# this is to say
	# time to fill is for the integral to equal to 1 (since everything equlas to 1)
	# \int_0^{time_to_fill}(d1+d2) dt = 1
	# therefore time_to_fill = 1/(d1+d2)
	# if t_max is greater than that, we fill, otherwise we get a fraction
	# note this does not take into accountcompeteing for the place in the
	# membrane, but rather takes that finding the spot is not the problem
	# the numerical difference is not that big, but the simulation, in this python version
	# fiftten time faster
	sumd = sum(delivery_rate)
	if sumd<0.000001: return [0.0, 0.0]
	time_to_fill = 1.0/sumd
	if max_steps>=time_to_fill:
		return [delivery_rate[0]/sumd, delivery_rate[1]/sumd]
	else:
		fraction = max_steps/time_to_fill
		return [delivery_rate[0]/sumd*fraction, delivery_rate[1]/sumd*fraction]

########################################
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
		free_space = (1.0 - total_fraction_occupied)
		# increase of each  population in the membrane  - this goes to zero as the free space does
		# this is different than the assumptions of filling_fractions_for_const_delivery_rates
		# it results in a little bit lower filling than filling_fractions_for_const_delivery_rates
		# free_space = 1 # check that the numbers are the same as in the filling_fractions_for_const_delivery_rates case
		delta = list(map(lambda a: a*free_space, available_pool))
		incorporated_fraction = [incorporated_fraction[i] + delta[i] for i in index_iterator]
		total_fraction_occupied = sum(incorporated_fraction)

		points.append(incorporated_fraction.copy())
		step +=1

	if verbose:
		print()
		print(f" %.3f   %.3f "%(delivery_rate[0], delivery_rate[1]))
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
			print("d1 relative %.3f"%(d1/(d1+d2)))
			print("d2 relative %.3f"%(d2/(d1+d2)))
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



########################################
def main():


	delivery_rate = [0.3, 0.5]
	for i in range(15):
		delivery_rate = [d/2 for d in delivery_rate]
		[a1, a2] = sim_core(delivery_rate, plot=False, verbose=False)
		[b1, b2] = filling_fractions_for_const_delivery_rates(delivery_rate)
		print(f" %.5f   %.5f        %.4f   %.4f      %.4f   %.4f "%(delivery_rate[0], delivery_rate[1], a1, b1, a2, b2))

	return

if __name__ == "__main__":
	main()