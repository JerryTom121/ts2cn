# code to analyse time series

# importing libraries
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, axes
from igraph import Graph
from igraph import ADJ_UNDIRECTED

# used to import cython code
#import pyximport
#pyximport.install()

 
# return the index of the peaks in the time series ts
def get_peaks(ts):
	"""
--------------------------------------------
Get the peaks of the time series
--------------------------------------------
ts:		List. The time series
--------------------------------------------
Return a list, each element represent index of the peak
--------------------------------------------
Usage example:

	"""
	peaks = []
	up = True
	down = False
	prev = ts[0]
	for i in range(1, len(ts)):
		if ts[i] < prev:
			down = True
		if ts[i] > prev:
			up = True
			down = False
		if up and down:
			peaks.append(i-1)
			down = False
			up = False
		prev = ts[i]
	return peaks

# returns the maximum correlation between the cycles
# TODO allow other way to compute the proximity
def proximity(x, y):
	"""
--------------------------------------------
Calculate the maximum correlation between cycles
--------------------------------------------
x:	List. Time series segments
y:	List. Time series segments
--------------------------------------------
Return a number. Calculate the maximu correlation by running the smaller segment over the bigger
--------------------------------------------
Usage example:

	"""
	from scipy.stats import pearsonr
	bigger = x
	smaller = y
	if len(x) < len(y):
		bigger = y
		smaller = x
	sm_size = len(smaller)
	bg_size = len(bigger)
	# in case the two cycles have the same length
	
	"""
	PROFILING
	from time import time
	t=time()
	if sm_size == bg_size:
		p= pearsonr(smaller, bigger)[0]
		print(('if', time()-t))
		return p
	p = max([pearsonr(smaller, bigger[i:i+sm_size])[0] for i in range(bg_size-sm_size)])
	print(('max',time()-t))
	return p
	"""
	
		
	if sm_size == bg_size:
		return pearsonr(smaller, bigger)[0]
	return max([pearsonr(smaller, bigger[i:i+sm_size])[0] for i in range(bg_size-sm_size)])
	

# TODO allow other way to "select" the cycles, currently cycle are peak to peak segments
# TODO improve the performance, it takes 70 seconds
def cycle_network(ts, cycle_ths=0.3):
	"""
--------------------------------------------
Convert a time series into a cycle network
--------------------------------------------
ts:			List. The time series
cycle_ths:	Number. The threshold of proximity between cycles
--------------------------------------------
Return a graph object, using the igraph representation
Currently the cycles are defined as the time series segments between peaks, and
the proximity measure between cycles is the pearson correlation

See article for details of the implementation
	Recurrence-Based Time Series Analysis by means of Complex Network Methods, 2011
	Donner, R.; Small, M.; Donges, J.; Marwan, N.; Zou, Y.; Xiang, R.; Kurths, J.
	International Journal of Bifurcation and Chaos, Vol. 21 No. 4 (2011)
--------------------------------------------
Usage example:

import numpy as np
import imp

from ts2cn.ts import cycle 

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

graph = cycle.cycle_network(ts, cycle_ths=0.3)
	"""
	from scipy.spatial.distance import squareform, pdist
	

	from time import time
	
	t = time()
	peaks = get_peaks(ts)
	cycles = np.array([ts[peaks[i-1]:peaks[i]] for i in range(1, len(peaks))])
	print(('get cycles', time()-t))


	t = time()
	# still has problem with types conversions
#	from ts2cn.ts.cycle_opt import cycle_opt
#	corr = cycle_opt.correlation(cycles, cycle_ths)


	# this is the one I will leave
#	corr = np.array([proximity(cycles[i], cycles[j]) for i in range(len(cycles)) for j in range(i+1, len(cycles))])

#	corr = pdist(cycles, metric=proximity)
#	corr = cycle_corr.corr(cycles)
#	corr = cycle_cython.corr(cycles)
	
	# this code has the best performance
	# Try this code as cython
	from scipy.stats import pearsonr
	#corr = []
	corr = [0] * int((len(cycles)*(len(cycles)-1))/2)
	l = 0
	for i in range(len(cycles)):
		for j in range(i+1, len(cycles)):
			bigger = cycles[i]
			smaller = cycles[j]
			bg_size = len(bigger)
			sm_size = len(smaller)
			if sm_size > bg_size:
				bigger = cycles[j]
				smaller = cycles[i]
			bg_size = len(bigger)
			sm_size = len(smaller)
			
			#connect = 0
			for k in range(bg_size - sm_size):
				r =  pearsonr(smaller, bigger[k:k+sm_size])[0]
				if r > cycle_ths:
					#connect = 1
					corr[l] = 1
					break
			l += 1
			#corr.append(connect)

	
	print(time()-t)
	
	adj_mat = squareform(corr)
	
#	adj_mat = squareform(corr > cycle_ths)
	return Graph.Adjacency(adj_mat.tolist(), mode=ADJ_UNDIRECTED)

