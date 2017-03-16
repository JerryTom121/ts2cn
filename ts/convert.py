# interface to the convertion methods

# importing libraries
import numpy as np

from igraph import Graph

# from phase_space
from ts2cn.ts.phase_space import reconstruct_ps, epsilon_recurrence_network, k_nearest_network, adaptative_nearest_network, correlation_network
# from transition_networks
from ts2cn.ts.transition import transition_network
from ts2cn.ts.cycle import cycle_network
# from visibility
from ts2cn.ts.visibility import visibility_network


# TODO allow for more parameters of the phase space reconstruction
# TODO avoid reconstructing the phase space several times
# interface for the complete analysis function
def ts_to_cn(ts, args, greater_component=True, min_nodes=5):
	"""
--------------------------------------------
Convert a time series into a complex network
--------------------------------------------
ts:					List. An array representing the time series
args:				Dict. The parameters of the conversion
					Possible keys:
						
						# cycle network parameter
						cycle_ths:	Number. The threshold of proximity between cycles

						
						# transition network parameters
						bins:		List. The numbers representing the bins separators
						prob_ths:	Number. The probability threshold, probability values above this value means 
									the points are connected
						
						# visibility graph parameter
						vibility:	The type of the visibility graph, "natural" or "horizontal"
						
						# related to the phase space reconstruction
						phase_space:		Array. Phase space in array form, it's is used to avoid the reconstruction several times. 
											See ts2cn.ts.phase_space.reconstruct_ps for more info.
						noise_perc:			Number. The percentage of noise allowed in finding the first minimum of the Average Mutual Information (AMI)
						false_nn_threshold:	Number. The percentagem of false nearest neighbors allowed during the unfolding process
						
						# epsilon recurrence network parameter
						epsilon:	Number or Dict. The distance threshold, distance values below this means the two 
									points should be connected. If a dict is passed instead, it should have the key 'percentile'
									with the corresponding percentile to compute.
						dist_vec:	Array. The distance matrix between points, if not passed this will be calculated. 
									Used to save computation if this is alread computed
						
						# k nearest neighbor network parameter
						k:		Number. The k nearest neighbors will be connected
						
						# adaptative nearest network parameter
						E0:		Number.  The E0 nearest neighbors will be connected according to specification
						
						# correlation network parameter
						rho_ths:	Number. The rho (greek letter used for correlation) threshold.
									Points with correlation above this value will be connected.
									This varies between -1 and +1
						
greater_component:	Boolean. True to extract the greater component of the graph
min_nodes:			Number. The minium number of nodes the phase space should have else a empty graph is returned
--------------------------------------------
Return a graph object, using igraph representation
All networks but the transition and k nearest neighbor networks are undirected

ATTENTION networks that require phase space reconstruction may used a lot of RAM memory

For more details of each network conversion look for help of the related network conversion method

See article for details of the implementations
	Recurrence-Based Time Series Analysis by means of Complex Network Methods, 2011
	Donner, R.; Small, M.; Donges, J.; Marwan, N.; Zou, Y.; Xiang, R.; Kurths, J.
	International Journal of Bifurcation and Chaos, Vol. 21 No. 4 (2011)

--------------------------------------------
Usage example:

import numpy as np
import imp

from ts2cn.ts import convert

filename='rossler.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

# cycle network with a proximity cycles threshold of 0.3
graph_cycle = convert.ts_to_cn(ts, args={'cycle_ths': 0.3}, greater_component=True)

# transition network with a probability threshold of 0.3
bins = np.linspace(start=min(ts)-1, stop=max(ts)+1, num=100)
graph_transition = convert.ts_to_cn(ts, args={'bins': bins, 'prob_ths':0.3}, greater_component=True)

# visibility graph, the natural visibility graph
graph_visibility = convert.ts_to_cn(ts, args={'visibility': 'natural'}, greater_component=True)

# all networks that require the phase space reconstruction use false_nn_threshold 0.2 and noise_perc 2

# epsilon recurrence network with epsilon being the 25-th percentile of the phace space distances
graph_epsilon = convert.ts_to_cn(ts, args={'epsilon': {'percentile': 25}, 'false_nn_threshold': 0.2, 'noise_perc': 2}, greater_component=True)

# k nearest network connecting the 5 closest points
graph_knn = convert.ts_to_cn(ts, args={'k': 5, 'false_nn_threshold': 0.2, 'noise_perc': 2}, greater_component=True)

# adaptative nearest network connecting the 5 closest points, with the adaptative approach
graph_ann = convert.ts_to_cn(ts, args={'E0': 5, 'false_nn_threshold': 0.2, 'noise_perc': 2}, greater_component=True)

# correlation network connecting point with a correlation above 0.2 (rho_ths)
graph_correlation = convert.ts_to_cn(ts, args={'rho_ths': 0.2, 'false_nn_threshold': 0.2, 'noise_perc': 2}, greater_component=True)
	"""
	# First the methods that do not require the phase space reconstruction
	
	directed = False
	if 'directed' in args:
		directed = args['directed']

	# cycle network
	if 'cycle_ths' in args:
		graph = cycle_network(ts, cycle_ths=args['cycle_ths'])
	# transition probability network
	elif 'bins' in args and 'prob_ths' in args:
		graph = transition_network(ts, bins=args['bins'], threshold=args['prob_ths']).simplify()
	# visibility network
	elif 'visibility' in args:
		graph = visibility_network(ts, vis_type=args['visibility'])
	# convert_method  and convert_args in the kwargs parameter

	else:
		# Methods that require the phase space reconstruction (recurrence)
		noise_perc = 2
		if 'noise_perc' in args:
			noise_perc = args['noise_perc']
		# avoid reconstruction of phase space several times
		if 'phase_space' in args:
			rc = args['phase_space']
		else:
			rc = reconstruct_ps(ts, false_nn_threshold=args['false_nn_threshold'], noise_perc=noise_perc)

		# avoid try to convert a phase space with less than min_nodes
		if len(rc) < min_nodes:
			return Graph()
		
		# epsilon-recurrent network
		if 'epsilon' in args:
			graph = epsilon_recurrence_network(rc, epsilon=args['epsilon'], dist_vec=None)
		# k-nearest-neighbor network
		elif 'k' in args:
			graph = k_nearest_network(rc, k=args['k'])
		# Adaptive-nearest-neighbor network
		elif 'E0' in args:
#			graph = adaptative_nearest_network(rc, E0=args['E0'])
			# without simplify some graphs apear to be not simple ones
			graph = adaptative_nearest_network(rc, E0=args['E0']).simplify()
		# correlation network
		elif 'rho_ths' in args:
			graph = correlation_network(rc, rho_ths=args['rho_ths'])

	# extract greater component
	if greater_component:
		return extract_greater_component(graph)
	return graph



# Extract the greater component of a given graph
def extract_greater_component(graph):
	from igraph import STRONG, WEAK
	cls = graph.components(mode=STRONG)
	idxs = [len(i) for i in cls]
	return cls.subgraph(idxs.index(max(idxs)))



