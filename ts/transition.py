# code to convert time series to transition network

# importing libraries
from igraph import *
from numpy import array, digitize, linspace, dtype
from scipy.sparse import coo_matrix, lil_matrix, csr_matrix

# TODO allow to pass the percentile of the probabilities like the epsilon recurrence network

# interface for the complete analysis function
def transition_network(ts, bins, threshold):
	"""
--------------------------------------------
Convert the time series into a transition network
--------------------------------------------
ts:			List. The time series
bins:		List. The numbers representing the bins separators
threshold:	Number. The probability threshold, probability values above this value means 
			the points are connected
--------------------------------------------
Return a graph object, using the igraph representation

ATTENTION beware of the first and last values of the bins, it should be a littee smalle and
a little bigger than the min and max value of time series


See article for details of the implementation
	Recurrence-Based Time Series Analysis by means of Complex Network Methods, 2011
	Donner, R.; Small, M.; Donges, J.; Marwan, N.; Zou, Y.; Xiang, R.; Kurths, J.
	International Journal of Bifurcation and Chaos, Vol. 21 No. 4 (2011)
--------------------------------------------
Usage example:

import numpy as np
from ts2cn.ts import transition as trs

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

bins = np.linspace(start=min(ts)-1, stop=max(ts)+1, num=100)

graph = trs.transition_network(ts, bins=bins, threshold=0.3)
	"""
	return Graph.Adjacency(
			matrix=(transition_prob(ts=ts, bins=bins) >= threshold).toarray().tolist(),
			mode=ADJ_DIRECTED)


# return the transition probability matrix
def transition_prob(ts, bins):
	"""
--------------------------------------------
Compute the transition probability matrix
--------------------------------------------
ts:		List. The time series
bins:	List. The numbers representing the bins separators
--------------------------------------------
Return a matrix with the probability, using scipy sparse matrix representation (CSR)
--------------------------------------------
Usage example:

import numpy as np
from ts2cn.ts import transition as trs

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

bins = np.linspace(start=min(ts), stop=max(ts), num=100)

probability_matrix = trs.transition_prob(ts, bins=bins)
	"""
	# create an empty sparse matrix 
	Cij = csr_matrix((len(bins)-1, len(bins)-1), dtype=dtype(float))
	
	# put each time series observation into his respective class (defined by the bins)
	# It starts at 1, so in order to maintain the indexes 1 is subtracted
	S = digitize(x=ts,bins=bins) -1

	# runs over S and count the transition between classes
	for i in range(len(S)-1):
		Si = S[i]
		Sj = S[i+1]
		Cij[Si,Sj] += 1
	
	return Cij/len(ts)

