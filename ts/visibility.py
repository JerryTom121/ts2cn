# code to convert time series to visibility graph

# importing libraries
from igraph import *
from numpy import array, arange

# Attempt to use cython to improve the performance
#import pyximport
#pyximport.install()
#from ts2cn.ts import visibility_mat as vis_mat

from time import time

# TODO - Check to garantee the correctness of the conversion
# TODO - Implement other visibility types (horizontal)
def visibility_network(ts, vis_type='natural'):
	"""
--------------------------------------------
Function to convert a time series into a visibility graph
--------------------------------------------
ts:			Time series in array format
vis_type:	The type of the visibility graph, "natural" or "horizontal"
--------------------------------------------
Return a visibility graph in the igraph representation of the time series
--------------------------------------------
Usage example:

import numpy as np
from ts2cn.ts import visibility as vis

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

graph = vis.visibility_network(ts, vis_type='natural')
	"""
#graph='natural', undirected=True):

	


	size = len(ts)
	
	# create an empty adjacency matrix
	adj_mat = array([[0]*size] * size, dtype='int')
	
	# cython attempt to improve performance
#	adj_mat = array( vis_mat.matrix(array(ts, dtype=float), int(size)) )
	#adj_mat = vis_mat.matrix( array(ts, dtype=float), int(size))

	# attempt to improve through pdist function
	# not a good choice, takes too much ram 
	"""
	from scipy.spatial.distance import pdist, squareform
	def vis_dist(ya, yb):
		pass
	cte = squareform(pdist([(i, 0) for i in ts]))
	# connecting the direct neighbors
	adj_mat[0][1] = 1
	adj_mat[size-1][size-2] = 1
	for i in range(1, size-1):
		adj_mat[i][i-1] = 1
		adj_mat[i][i+1] = 1
	
	yc_tc = [(ts[(i+1):(j-1)], arange(start=i+1, stop=j-1)) for i in range(size-1) for j in range(i+2, size)]

	for i in range(size-1):
		for j in range(i+2, size):
#			yc = ts[(i+1):(j-1)]
#			tc = arange(start=i+1, stop=j-1)
			# do not connect i and j if the condition is satisfied
			#if (yc + tc*cte < ts[j] + cte*j).all():
			if (yc_tc[i,j][0] <  ts[j] + cte[i,j]*(j - tc[i,j][1])).all():
				continue
			adj_mat[i][j], adj_mat[j][i] = 1, 1


	t = time()
	total = time()
	"""

	t = time()
	total = time()
	# connecting the direct neighbors
	adj_mat[0][1] = 1
	adj_mat[size-1][size-2] = 1
	for i in range(1, size-1):
		adj_mat[i][i-1] = 1
		adj_mat[i][i+1] = 1
	
	print('first '+ str(time() - t))
	t = time()
#	""
	for i in range(size-1):
		for j in range(i+2, size):
			cte = (ts[i]-ts[j])/(j-i)
			yc = ts[(i+1):(j-1)]
			tc = arange(start=i+1, stop=j-1)
			# do not connect i and j if the condition is satisfied
			#if (yc + tc*cte < ts[j] + cte*j).all():
			if (yc <  ts[j] + cte*(j - tc)).all():
				continue
			adj_mat[i][j], adj_mat[j][i] = 1, 1
#	"""
#	a = [(i,j) for i in range(size-1) for j in range(i+2, size) if not (ts[(i+1):(j-1)] + arange(start=i+1, stop=j-1)*(ts[i]-ts[j])/(j-i) < ts[j] + (ts[i]-ts[j])/(j-i)*j ).all()]
			
	
	# TODO - implement other types of visibility graphs (horizontal and directed)
	
	"""
	# iterate over the elements of the time series
	for i in range(size-1):
		
		# set the "maximum" value for elements on the right of the current value i
		k = i+1
		
		# it starts in i+2 because every direct neighbor is connected
		# verify which elements on the right of the current element i,
		# the j elements, should be connected with i
		for j in range(i+2, size):
			
			# TODO - test another algebraic inequality
			
			# if the inequality is true connect i and j
			if ts[k] < ts[j] + (ts[i]-ts[j])*(j-k)/(j-i):
				adj_mat[i][j] = 1
				adj_mat[j][i] = 1
			else:
				k = j
#	"""
	
	
	print('main '+ str(time() - t))
	t = time()

	print('total '+ str(time() - total))
	
#	igraph.plot(igraph.Graph.Adjacency([[1,1,1],[1,1,0],[0,0,1]]))
	return Graph.Adjacency(adj_mat.tolist(), mode=ADJ_UNDIRECTED)

