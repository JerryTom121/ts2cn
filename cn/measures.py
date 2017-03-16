# code to calculate the desired measures from the complex network

# importing libraries
from igraph import *
from numpy import array, nan

#####
# Rule of thumb of this implementation, if there is a igraph method
# that do the job it will be used
####

def dist_n_moment(dist=[], n=2):
	"""
--------------------------------------------
Calculates the n-th distribution moment
--------------------------------------------
dist:	The distribuition in igraph's Histogram form
n:	The number of the desired moment
--------------------------------------------
Returns a number
--------------------------------------------
	"""
	# i[2]*(i[1]-i[0]) gives the probability of the given bin
	# i[1]-i[0] is the bin's width, and i[2] the bin's height
	# i[0] is the reference value of the bin and the value
	# to be raised in order to calculate the moment
	return sum([i[2]*(i[1]-i[0])*i[0]**n for i in dist])	


def degree(graph, vertices=None):
	"""
--------------------------------------------
Degree 
--------------------------------------------
graph:	An igraph object 
vertices:	A list, the vertices to return the degree
--------------------------------------------
Return a list
--------------------------------------------
	"""
	return graph.degree(vertices=vertices)


def degree_distribution(graph, freq=True, *args, **kwargs):
	"""
--------------------------------------------
Degree distribution
--------------------------------------------
graph:	An igraph object 
freq:		Boolean, true to return the frequency 
			instead of the counting in the bins
--------------------------------------------
Returns a dictionary with keys 'mean', 'var', 'sd', 'n' and 'bins'
The key 'bins' hold a list with the bins, it does not necessarily 
starts at the value 1, the first bin is the one with the first
counting
--------------------------------------------
	"""
	dist = graph.degree_distribution()
	
	if freq:
		dist.bins = [(i[0], i[1], i[2]/dist.n) for i in dist.bins()]
	else:
		dist.bins = [i for i in dist.bins()]
	
	return {
		'mean': dist.mean,
		'var': dist.var,
		'sd': dist.sd,
		'n': dist.n,
		'bins': dist.bins
		}
def degree_entropy(dist=[], normalize=False):
	"""
--------------------------------------------
Calculate the degree entropy
--------------------------------------------
dist:			A list representing the distribution, bins, 
				from which the entropy will be calculated
normalize:    Boolean, if normalize the entropy
--------------------------------------------
Returns a number
--------------------------------------------
	"""
	from math import log
	entropy = 0

	for k in range(len(dist)):
		if dist[k][2]:
			prob_k = (dist[k][1] - dist[k][0]) * dist[k][2]
			entropy -= prob_k  * log(prob_k)

	if normalize:
		# normalization means to divide the entropy by the 
		# total number of bins, i.e, the greatest degree found
		return entropy / log(dist[len(dist)-1][0])
	return entropy


def clustering_coef(graph, vertices=None):
	"""
--------------------------------------------
Calculate the clustering coefficient
--------------------------------------------
graph:	An igraph object
vertices:	List or number of the desired vertices to calculate
			the clustering coefficient
--------------------------------------------
Returns a list with the clustering coefficient of each vertex
or a number if only one vertex is asked
--------------------------------------------
	"""
	return graph.transitivity_local_undirected(vertices=vertices)

def transitivity(graph):
	"""
--------------------------------------------
Calculate the graph transitivity (global clustering coefficient)
-------------------------------------------
graph:	An igraph object
--------------------------------------------
Returns a number
--------------------------------------------
	"""
	return graph.transitivity_undirected()


def efficiency(graph):
	"""
--------------------------------------------
Calculates the graph Efficiency
--------------------------------------------
graph:	An igraph object
--------------------------------------------
Returns a number
According to Efficient Behavior of Small-World Networks
Vito Latora and Massimo Marchiori
Phys. Rev. Lett. 87, 198701 – Published 17 October 2001
--------------------------------------------
	"""
	n = len(graph.vs())
	return sum([sum([ 1/j for j in i if j]) for i in graph.shortest_paths()]) / (n*(n-1))

def complexity(dist=[]):
	"""
--------------------------------------------
Calculate the graph Complexity
--------------------------------------------
dist:	A list with the degree distribution, in histogram form
--------------------------------------------
Returns a number
It measures how connected is the graph
--------------------------------------------
	"""
	mean = dist_n_moment(dist=dist, n=1)
	return dist_n_moment(dist=dist, n=2) / mean if mean != 0 else nan
#	return dist_n_moment(dist=dist, n=2) / mean 

def knn(graph, vertices=None):
	"""
--------------------------------------------
Calculate the graph average degree of the neighbors for each vertex
--------------------------------------------
graph:	An igraph object
vertices:	A list of vertices.
--------------------------------------------
Returns a list 
The average degree of neighbors as a function of vertex degree. 
The zeroth element of this list corresponds to vertices of degree 1.
Not-a-number (nan) can be returned if the value there is not defined
--------------------------------------------
	"""
#	plot(graph)
#	print(vertices)
#	print(graph.simplify().knn()[0])
#	print(graph.simplify().knn()[1])
	return graph.knn(vids=vertices)[1]


def assortativity(graph, directed=False):
	"""
--------------------------------------------
Calculate the graph assortativity
-------------------------------------------
graph:	An igraph object
directed:	Boolean, tells if the graph is directed or not
--------------------------------------------
Returns a number 
--------------------------------------------
	"""
	return graph.assortativity_degree(directed=directed)


def central_point_dominance(betweenness=[]):
	"""
--------------------------------------------
Calculate the graph's central point dominance
--------------------------------------------
betweenness:	A list with the betweenness of each vertex
--------------------------------------------
Returns a number 
--------------------------------------------
	"""
	bt_max = max(betweenness)
	return sum([bt_max - i for i in betweenness]) / len(betweenness)


def betweenness(graph, vertices=None, directed=False):
	"""
--------------------------------------------
Calculate the betweenness centrality
--------------------------------------------
graph:	An igraph object
vertices:	A list of vertices
directed:	Boolean, tells if the graph is directed or not
--------------------------------------------
Returns a list with the betweenness of each vertex
--------------------------------------------
	"""
	return graph.betweenness(vertices=vertices, directed=directed)


def eigenvector(graph, directed=False, normalize=True):
	"""
--------------------------------------------
Calculate the eigenvector centrality
--------------------------------------------
graph:		An igraph object
directed:		Boolean, tells if the graph is directed or not
normalize:	Boolean, if the result should be normalized
--------------------------------------------
Returns a list with the eigenvector centrality of each vertex
If normalized the largest one will be 1
--------------------------------------------
	"""
	return graph.eigenvector_centrality(directed=directed, scale=normalize, 
		weights=None, return_eigenvalue=False)

def pagerank(graph, vertices=None, directed=False, damping=0.85):
	"""
--------------------------------------------
Calculate the PageRank centrality
--------------------------------------------
graph:	An igraph object
vertices:	A list of vertices
directed:	Boolean, tells if the graph is directed or not
damping:	A number, the damping factor	
--------------------------------------------
Returns a list with the PageRank centrality of each vertex
--------------------------------------------
	"""
	return graph.pagerank(vertices=vertices, directed=directed, damping=damping, 
		weights=None, implementation='prpack', niter=1000, eps=0.001)

def closeness(graph, vertices=None, normalize=True, mode='ALL'):
	"""
--------------------------------------------
Calculate the closeness centrality
--------------------------------------------
graph:		An igraph object
vertices:		A list of vertices
normalize:	Boolean, if the result should be normalized
mode:			String, the mode of the coreness, 'ALL', 'IN' or 'OUT'
--------------------------------------------
Returns a list with the closeness centrality of each vertex
--------------------------------------------
	"""
	return graph.closeness(vertices=vertices, normalized=normalize, 
		weights=None, cutoff=None)

def coreness(graph, mode='ALL'):
	"""
--------------------------------------------
Calculate the coreness (K-core) centrality
--------------------------------------------
graph:	An igraph object
mode:		String, the mode of the coreness, 'ALL', 'IN' or 'OUT'
--------------------------------------------
Returns a list with the coreness centrality of each vertex
--------------------------------------------
	"""
	return graph.coreness(mode=mode)





# TODO calculate the power law fit
# TODO maybe it should be the case of see a goodness of fit for the normal for the ditrbutions calculated
def extract_metrics(graph, aggregate=True, min_nodes=5):
	"""
--------------------------------------------
Extract the basic metrics of the complex network
---------------------------------------------
graph:      An igraph object
aggregate:  Bollean, if the histograms should be summarized
            returning only the mean and standard deviation
min_nodes:	Number. The minimum number of nodes the graph needs
			in order to compute the metrics. If the graph has less
			than this number, then nan (not a number) is return
			for all metrics
--------------------------------------------
Return a dictionary with the following fields:
degree      degree - list
dg_entropy  degree entropy - number
clus_coef   clustering coefficient - list
trans       graph's transitivity - number
eff         efficiency - number
compl       complexity - number
knn         mean k-nearest-neighboor as function of the degree - list
assort      assortativity - number
betweenness betweenness centrality - list
cpd         central point dominance - number
eigenvector eigenvector centrality - list
pagerank    PageRank centrality - list
closeness   Closeness centrality - list
coreness    Coreness centrality - list
--------------------------------------------
	"""
	from math import isnan
	from numpy import average, std, nan, isnan, histogram, sqrt
	
	print(('NODES ', len(graph.vs)))
	
	# Asks for a minimum number of nodes to avoid
	# problem when computing some metrics (knn for instance)
	if len(graph.vs) < min_nodes:
		nan_return = {'mean': nan, 'var':nan}
		return {
			'degree': nan_return,
			'dg_entropy': nan,
			'clus_coef': nan_return,
			'trans': nan,
		#       'eff': nan,
			'compl': nan,
			'knn': nan_return,
			'assort': nan,
			'betweenness': nan_return,
			'cpd': nan,
			'eigenvector': nan_return,
			'pagerank': nan_return,
			'closeness': nan_return,
			'coreness': nan_return
		}

	# general metrics
	dg = degree(graph)

	dist = Histogram()
	dist << dg
	dg_dist = [(i[0], i[1], i[2]/dist.n) for i in dist.bins()]
	#   dist = None

	#dist = histogram(dg, bins=max(dg))
	#dg_dist = []

	dg_entropy = degree_entropy(dg_dist)
	clus_coef = clustering_coef(graph)
	trans = transitivity(graph)
	#eff = efficiency(graph)
	compl = complexity(dg_dist)
	k_nn = knn(graph)
	assort = assortativity(graph)

	# centralities
	bt_ness = betweenness(graph)
	cpd = central_point_dominance(bt_ness)
	egn_vector = eigenvector(graph)
	pg_rank = pagerank(graph)
	cl_ness = closeness(graph)
	core_ness = coreness(graph)

	# there is a method on igraph to fit the degree distribution to a power law
	# igraph.power_law_fit and it seems to be fast, good
	# just need to check its correctness.
	# Need to execute with some previous code and see if the results check

	# after the power law fit aggregate histograms if required
	if aggregate:
		"""
		objs_to_aggregate = [degree, clus_coef, knn, betweenness, 
			eigenvector, pagerank, closeness, coreness]
		for i in range(len(objs_to_aggregate)):
			dist = Histogram()
			# remove NaN elements, so far only seen in the knn array
			dist << [j for j in objs_to_aggregate[i] if not isnan(j)]
			objs_to_aggregate[i] = {'mean': dist.mean, 'var': dist.var}
		"""
		dg = {'mean': dist.mean, 'var': sqrt(dist.var)}

		#dist = Histogram()
		#dist << clus_coef
		#clus_coef = {'mean': dist.mean, 'var': dist.var}
		clus_coef = {'mean': average(clus_coef), 'var': std(clus_coef)}

		#dist = Histogram()
		# remove NaN elements, so far only seen in the k_nn array
		#dist << [j for j in k_nn if not isnan(j)]
		#k_nn = {'mean': dist.mean, 'var': dist.var}Vyp
		not_nan = [j for j in k_nn if not isnan(j)]
		k_nn = {'mean': average(not_nan), 'var': std(not_nan)}

		#dist = Histogram()
		#dist << bt_ness
		#bt_ness = {'mean': dist.mean, 'var': dist.var}
		bt_ness = {'mean': average(bt_ness), 'var': std(bt_ness)}

		#dist = Histogram()
		#dist << egn_vector
		#egn_vector = {'mean': dist.mean, 'var': dist.var}
		egn_vector = {'mean': average(egn_vector), 'var': std(egn_vector)}

		#dist = Histogram()
		#dist << pg_rank
		#pg_rank = {'mean': dist.mean, 'var': dist.var}
		pg_rank = {'mean': average(pg_rank), 'var': std(pg_rank)}

		#dist = Histogram()
		#dist << cl_ness
		#cl_ness = {'mean': dist.mean, 'var': dist.var}
		cl_ness = {'mean': average(cl_ness), 'var': std(cl_ness)}

		#dist = Histogram()
		#dist << core_ness
		#core_ness = {'mean': dist.mean, 'var': dist.var}
		core_ness = {'mean': average(core_ness), 'var': std(core_ness)}



	return {
		'degree': dg,
		'dg_entropy': dg_entropy,
		'clus_coef': clus_coef,
		'trans': trans,
	#       'eff': eff,
		'compl': compl,
		'knn': k_nn,
		'assort': assort,
		'betweenness': bt_ness,
		'cpd': cpd,
		'eigenvector': egn_vector,
		'pagerank': pg_rank,
		'closeness': cl_ness,
		'coreness': core_ness
}


