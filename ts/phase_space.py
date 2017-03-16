# code to convert time series to networks related to phase space

# importing libraries
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, axes
from igraph import Graph
from igraph import ADJ_UNDIRECTED


def reconstruct_ps(ts, max_dim=20, dims_step=5, false_nn_threshold=0.01, noise_perc=2):
	"""
--------------------------------------------
Reconstruct the time series's phase space
--------------------------------------------
ts:					List. The time series
max_dim:			Number. The maximum dimension that the reconstruction process will try
dims_step:			Number. The amount of dimensions that the reoncstruction process will try on each round. This is used to speed up the reconstruction
false_nn_threshold:	Number. The percentagem of false nearest neighbors allowed during the unfolding process
noise_perc:			Number. The percentage of noise allowed in finding the first minimum of the Average Mutual Information (AMI)
--------------------------------------------
Return an array, where each element represents a point in the phase space
--------------------------------------------
Usage example:

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import imp

from ts2cn.ts import phase_space as phs


filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]
plt.plot(ts)
plt.show()
plt.plot(phs.lagged_ami(ts, min_lag=0, max_lag=len(ts)/2)[1])
plt.show()


rc = phs.reconstruct_ps(ts, max_dim=20, dims_step=5, false_nn_threshold=0.2, noise_perc=2)

imp.reload(plt)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rc[:, 0], rc[:, 1], rc[:, 2])
plt.show()
	"""
	
	# compute the lag given the percentage of noise in the mutual information
	lag = compute_lag(ts, lag_step=1, n_bins=10, noise_perc=2)
	print(('lag', lag))

	# calculate the embedding dimensions
	# TODO receive a parameter as the threshold to determine the optimum embedding dimension
	#max_dim = 100
	#dims_step = 10
	dims = []
	#false_nn_threshold = 0.01
	for min_dim in range(1, max_dim, dims_step):
		dims = global_false_nearest_neighbors(ts, lag=lag, min_dims=min_dim, max_dims=min_dim+dims_step-1)
		print(('dims', dims))
		if min(dims[1]) <= false_nn_threshold:
			dims = list(dims[1] < false_nn_threshold).index(True) + min_dim
			break
	embedding_dim = dims
	
	print(('embedding_dim', embedding_dim))
	# TODO return the other calculated values 

	return reconstruct(ts, lag, n_dims=embedding_dim)

"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import imp

from ts2cn.ts import phase_space as phs



filename='lorenz.dat'; file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r'); ts = file.read().split(); ts=[float(i) for i in ts]; plt.plot(ts); plt.show();plt.plot(phs.lagged_ami(ts, min_lag=0, max_lag=len(ts)/2)[1]); plt.show()


rc = phs.reconstruct_ps(ts, false_nn_threshold=0.2)

imp.reload(plt)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(rc[:, 0], rc[:, 1], rc[:, 2]); plt.show()





import numpy as np
import matplotlib as mpl
mpl.use('pdf')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import imp

from ts2cn.ts import phase_space as phs

filename='exptchao.dat'; file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r'); ts = file.read().split(); ts=[float(i) for i in ts]; plt.plot(ts); plt.show();

rc = phs.reconstruct_ps(ts, false_nn_threshold=0.2)

imp.reload(plt)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rc[:, 0], rc[:, 1], rc[:, 2]); plt.savefig('/tmp/figs/'+filename+'_st.png')

import time
from scipy.spatial.distance import pdist, squareform

imp.reload(plt); plt.close();
ts=time.time(); dist=pdist(rc); print(time.time()-ts); ts=time.time(); mat=squareform(dist); print(time.time()-ts); ts=time.time(); plt.imshow(mat); plt.savefig('/tmp/figs/'+filename+'_heat.png'); print(time.time()-ts);






from pilot.pilot import teste1

tss = teste1()[0]
ts = tss[1]['S'][0]
plt.plot(ts); plt.show()

rc = phs.reconstruct_ps(ts, false_nn_threshold=0.1)

imp.reload(plt)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rc[:, 0], rc[:, 1], rc[:, 2]); plt.show()


"""



# TODO incorporate somehow a measure to get a good value for the thresholds
# like analyse the histogram and get the value related to the first quartile for instance
"""
bins = 10
hh = np.histogram(data, density=True, bins=bins)
qt = 0.25
qt /= hh[1][1] - hh[1][0]
density_sum = 0
ths = hh[0][0]
for i in enumerate(hh[0]):
	density_sum += i[1]
	if density_sum > qt:
		ths = hh[1][ i[0] ]
		break


# MAIS FACIL, retorna o qt-percentile dos dados
qt = 25
ths = np.percentile(a=data, q=qt)
"""

# TODO check the use of sparse matrix in order to save RAM
def epsilon_recurrence_network(phase_space, epsilon=0.1, dist_vec=None):
	"""
--------------------------------------------
Convert a phase space into a epsilon recurrence network, an undirected one
--------------------------------------------
phase_space:	Array. The phase space representation in numpy array format
epsilon:		Number or Dict. The distance threshold, distance values below this means the two 
				points should be connected. If a dict is passed instead, it should have the key 'percentile'
				with the corresponding percentile to compute.
dist_vec:		Array. The distance matrix between points, if not passed this will be calculated. 
				It's used to save computation if this is alread computed
--------------------------------------------
Return a graph object, using igraph representation
--------------------------------------------
Usage example:

import numpy as np
import imp

from ts2cn.ts import phase_space as phs

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

rc = phs.reconstruct_ps(ts, max_dim=20, dims_step=5, false_nn_threshold=0.2, noise_perc=2)

graph = phs.epsilon_recurrence_network(rc, {'percentile': 25})
	"""
	from scipy.spatial.distance import pdist, squareform
	from igraph import Graph
	from igraph import ADJ_UNDIRECTED
	
	# TODO allow other distance metrics
	if not dist_vec:
		dist_vec = pdist(X=phase_space, metric='euclidean')
	
	# if set to get the epsilon automatic, get the value corresponding to 
	# the passed percentile of the distances
	if type(epsilon) == type({}):
		epsilon = np.percentile(dist_vec, epsilon['percentile'])
	
	adj_mat = squareform( dist_vec < epsilon )
	return Graph.Adjacency(adj_mat.tolist(), mode=ADJ_UNDIRECTED)
	



# TODO check the use of sparse matrix in order to save RAM
def k_nearest_network(phase_space, k=5):
	"""
--------------------------------------------
Convert a phase space into a k nearest neighbor network, a directed one
--------------------------------------------
phase_space:	Array. The phase space representation in numpy array format
k:				Number. The k nearest neighbors will be connected
--------------------------------------------
Return a graph object, using igraph representation
--------------------------------------------
Usage example:

import numpy as np
import imp

from ts2cn.ts import phase_space as phs

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

rc = phs.reconstruct_ps(ts, max_dim=20, dims_step=5, false_nn_threshold=0.2, noise_perc=2)

graph = phs.k_nearest_network(rc, k=5)

	"""
	from scipy.spatial.distance import pdist, squareform
	from igraph import Graph
	from igraph import ADJ_UNDIRECTED, ADJ_DIRECTED
	from sklearn.neighbors import NearestNeighbors

	# TODO allow other algorithms
	# it's passed k+1 because each node is considered the nearest neighboor of itself
	nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='kd_tree').fit(phase_space)
	
	adj_mat = nbrs.kneighbors_graph(phase_space, mode='connectivity').toarray()
	diag = range(len(adj_mat))
	adj_mat[diag, diag] = 0
	return Graph.Adjacency(adj_mat.tolist(), mode=ADJ_DIRECTED)
	
"""
PROBLEMA é descobrir qual o k adequado
No link sobre o knn do scipy ele fala de um possível problema da kd_tree para alta dimensões,
mas no caso as dimensões não passa de 7 talvez 8 então creio não ser um problema
"""


# TODO check the use of sparse matrix in order to save RAM
def adaptative_nearest_network(phase_space, E0=5):
	"""
--------------------------------------------
Convert a phase space into a adaptative nearest neighbor network, an undirected one
--------------------------------------------
phase_space:	Array. The phase space representation in numpy array format
E0:				Number.  The E0 nearest neighbors will be connected according to specification
--------------------------------------------
Return a graph object, using igraph representation
See article for details of the implementation
	Recurrence-Based Time Series Analysis by means of Complex Network Methods, 2011
	Donner, R.; Small, M.; Donges, J.; Marwan, N.; Zou, Y.; Xiang, R.; Kurths, J.
	International Journal of Bifurcation and Chaos, Vol. 21 No. 4 (2011)
--------------------------------------------
Usage example:

import numpy as np
import imp

from ts2cn.ts import phase_space as phs

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

rc = phs.reconstruct_ps(ts, max_dim=20, dims_step=5, false_nn_threshold=0.2, noise_perc=2)

graph = phs.adaptative_nearest_network(rc, E0=5)

	"""
	from scipy.spatial.distance import pdist, squareform
	from igraph import Graph
	from igraph import ADJ_UNDIRECTED, ADJ_DIRECTED
	from sklearn.neighbors import NearestNeighbors

	# TODO allow other algorithms
	# it's passed E0+1 because each node is considered the nearest neighboor of itself
	nbrs = NearestNeighbors(n_neighbors=E0+1, algorithm='kd_tree').fit(phase_space)
	
	adj_mat = nbrs.kneighbors_graph(phase_space, mode='connectivity')
	
	diag = range(adj_mat.shape[0])
	adj_mat[diag, diag] = 0
	
	# adaptation of edges
	for i in diag:
		k_offset = adj_mat[0:i, i].sum() +1
		# discard the first neighbor, itself
		new_neighbors = nbrs.kneighbors(X=[phase_space[i]], n_neighbors=k_offset, return_distance=False)[0][1:]
		adj_mat[i, new_neighbors] = 1
	
	return Graph.Adjacency(adj_mat.toarray().tolist(), mode=ADJ_UNDIRECTED)
	
"""
PROBLEMA é descobrir qual o k adequado
No link sobre o knn do scipy ele fala de um possível problema da kd_tree para alta dimensões,
mas no caso as dimensões não passa de 7 talvez 8 então creio não ser um problema
"""


# TODO improve the performance, it takes 1250 seconds
# maybe call the R code to perform the correlation, it is faster, it took in python
# 70 seconds for a 1000 elements array and 50 in R for the same array
# It can be done by saving the array in a file and then call an R code to read and process it
# TODO check the use of sparse matrix in order to save RAM
def correlation_network(phase_space, rho_ths=0.2):
	"""
--------------------------------------------
Convert a phase space into a correlation network, an undirected one
--------------------------------------------
phase_space:	Array. The phase space representation in numpy array format
rho_ths:		Number. The rho (greek letter used for correlation) threshold.
				Points with correlation above this value will be connected.
				This varies between -1 and +1
--------------------------------------------
Return a graph object, using igraph representation
See article for details of the implementation
	Recurrence-Based Time Series Analysis by means of Complex Network Methods, 2011
	Donner, R.; Small, M.; Donges, J.; Marwan, N.; Zou, Y.; Xiang, R.; Kurths, J.
	International Journal of Bifurcation and Chaos, Vol. 21 No. 4 (2011)
--------------------------------------------
Usage example:

import numpy as np
import imp

from ts2cn.ts import phase_space as phs

filename='lorenz.dat'
file = open('ts2cn/thirdy_parties/minfo/data/'+filename, 'r')
ts = file.read().split()
ts = [float(i) for i in ts]

rc = phs.reconstruct_ps(ts, max_dim=20, dims_step=5, false_nn_threshold=0.2, noise_perc=2)

graph = phs.correlation_network(rc, rho_ths=0.2)

	"""
	from scipy.stats import pearsonr as r
	from scipy.spatial.distance import squareform, pdist
	
#	from time import time
#	t=time()
	
	# Garantee that this approach is equivalent to the previous one
	
	# the definition of a distance is 1 - corr(x,y), so in order to get the correlation
	# again it's needed to subtract one and multiple by -1
	corr = (pdist(phase_space, metric='correlation' ) -1)*-1
#	corr = pdist(phase_space, metric=lambda x, y: r(x,y)[0] )
	#corr = [r(phase_space[i], phase_space[j])[0] for i in range(len(phase_space)) 
	#												for j in range(i+1, len(phase_space))]
#	corr = [True if r(phase_space[i], phase_space[j])[0] >= rho_ths else False for i in range(len(phase_space)) 
#													for j in range(i+1, len(phase_space))]
#	print(('corr', time()-t))
#	t=time()
	corr_ths = np.array(corr) >= rho_ths
#	corr_ths = corr
#	print(('corr_ths', time()-t))
#	t=time()
	adj_mat = squareform( corr_ths )
#	print(('square', time()-t))
	
#	adj_mat = squareform(np.array([r(phase_space[i], phase_space[j])[0] for i in range(len(phase_space)) 
#													for j in range(i, len(phase_space))]) >= rho_ths )
	
	'''
	# code calling R to compute the pearson correlation
	import os, rpy2.robjects as ro
	# first save phase_space vector to a temporary file
	tmp_file = '/tmp/tt'
	sep = ';'
	size = 100
	phase_space[:size].tofile(tmp_file, sep=sep)
	r_script = """
	mat = matrix(as.numeric(read.table(file="%s", sep="%s")), nrow=%d, ncol=%d, byrow=TRUE)
	size = %d
	for (i in 1:(size-1)) for (j in (i+1):size) print(cor(mat[i,], mat[j,]))
	""" % (tmp_file, sep, size, len(phase_space[0]), size)
	ro.r(r_script)
	'''

	return Graph.Adjacency(adj_mat.tolist(), mode=ADJ_UNDIRECTED)
	









from operator import sub

import numpy as np
from sklearn import metrics
from sklearn.neighbors import NearestNeighbors
from toolz import curry


def global_false_nearest_neighbors(x, lag, min_dims=1, max_dims=10, **cutoffs):
    """
    Across a range of embedding dimensions $d$, embeds $x(t)$ with lag $\tau$, finds all nearest neighbors,
    and computes the percentage of neighbors that that remain neighbors when an additional dimension is unfolded.
    See [1] for more information.

    Parameters
    ----------
    x : array-like
        Original signal $x(t).
    lag : int
        Time lag $\tau$ in units of the sampling time $h$ of $x(t)$.
    min_dims : int, optional
        The smallest embedding dimension $d$ to test.
    max_dims : int, optional
        The largest embedding dimension $d$ to test.
    relative_distance_cutoff : float, optional
        The cutoff for determining neighborliness,
        in distance increase relative to the original distance between neighboring points.
        The default, 15, is suggested in [1] (p. 41).
    relative_radius_cutoff : float, optional
        The cutoff for determining neighborliness,
        in distance increase relative to the radius of the attractor.
        The default, 2, is suggested in [1] (p. 42).

    Returns
    -------
    dims : ndarray
        The tested dimensions $d$.
    gfnn : ndarray
        The percentage of nearest neighbors that are false neighbors at each dimension.

    See Also
    --------
    reconstruct

    References
    ----------
    [1] Arbanel, H. D. (1996). *Analysis of Observed Chaotic Data* (pp. 40-43). New York: Springer.

    """
    x = _vector(x)

    dimensions = np.arange(min_dims, max_dims + 1)
    false_neighbor_pcts = np.array([_gfnn(x, lag, n_dims, **cutoffs) for n_dims in dimensions if len(x)/n_dims/lag -2 > 0])
    return dimensions, false_neighbor_pcts


def _gfnn(x, lag, n_dims, **cutoffs):
    # Global false nearest neighbors at a particular dimension.
    # Returns percent of all nearest neighbors that are still neighbors when the next dimension is unfolded.
    # Neighbors that can't be embedded due to lack of data are not counted in the denominator.
    offset = lag*n_dims
    is_true_neighbor = _is_true_neighbor(x, _radius(x), offset)
    return np.mean([
        not is_true_neighbor(indices, distance, **cutoffs)
        for indices, distance in _nearest_neighbors(reconstruct(x, lag, n_dims))
        if (indices + offset < x.size).all()
    ])


def _radius(x):
    # Per Arbanel (p. 42):
    # "the nominal 'radius' of the attractor defined as the RMS value of the data about its mean."
    return np.sqrt(((x - x.mean())**2).mean())


@curry
def _is_true_neighbor(
        x, attractor_radius, offset, indices, distance,
        relative_distance_cutoff=15,
        relative_radius_cutoff=2
):
    distance_increase = np.abs(sub(*x[indices + offset]))
    return (distance_increase / distance < relative_distance_cutoff and
            distance_increase / attractor_radius < relative_radius_cutoff)


def _nearest_neighbors(y):
    """
    Wrapper for sklearn.neighbors.NearestNeighbors.
    Yields the indices of the neighboring points, and the distance between them.

    """
    distances, indices = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(y).kneighbors(y)
    for distance, index in zip(distances, indices):
        yield index, distance[1]


def reconstruct(x, lag, n_dims):
    """Phase-space reconstruction.

    Given a signal $x(t)$, dimensionality $d$, and lag $\tau$, return the reconstructed signal
    \[
        \mathbf{y}(t) = [x(t), x(t + \tau), \ldots, x(t + (d - 1)\tau)].
    \]

    Parameters
    ----------
    x : array-like
        Original signal $x(t)$.
    lag : int
        Time lag $\tau$ in units of the sampling time $h$ of $x(t)$.
    n_dims : int
        Embedding dimension $d$.

    Returns
    -------
    ndarray
        $\mathbf{y}(t)$ as an array with $d$ columns.

    """
    #x = _vector(x)

    #if lag * (n_dims - 1) >= x.shape[0] // 2:
    #    raise ValueError('longest lag cannot be longer than half the length of x(t)')

    x = np.squeeze(x)

    lags = lag * np.arange(n_dims)
    return np.vstack(x[lag:lag - lags[-1] or None] for lag in lags).transpose()




def compute_lag(ts, lag_step=1, n_bins=10, noise_perc=2):
	# calculate the mutual information
	mut_inf = lagged_ami(ts, min_lag=0, max_lag=int(len(ts)/2), lag_step=lag_step, n_bins=n_bins)
	
	# get the first minimum and set it as the lag 
	
	# TODO maybe improve in order not to get a minimum due to noise
	# this happens with the henon data, there is no minimun consider min_mut = 1
	min_mut = mut_inf[1][0]
	
	noise_level = (max(mut_inf[1]) - min(mut_inf[1])) /100 * noise_perc

	for mi in mut_inf[1][1:]:
		if mi < min_mut :
			min_mut = mi
		elif abs(mi - min_mut) > noise_level:
			break
	lag = list(mut_inf[1]).index(min_mut)
	if lag == len(mut_inf[1])-1 or list(mut_inf[1]).index(mi) == len(mut_inf[1])-1:
		lag = 1
	
	return lag



def ami(x, y=None, n_bins=10):
    """Calculate the average mutual information between $x(t)$ and $y(t)$.

    Parameters
    ----------
    x : array-like
    y : array-like, optional
        $x(t)$ and $y(t)$.
        If only `x` is passed, it must have two columns;
        the first column defines $x(t)$ and the second $y(t)$.
    n_bins : int
        The number of bins to use when computing the joint histogram.

    Returns
    -------
    scalar
        Average mutual information between $x(t)$ and $y(t)$, in nats (natural log equivalent of bits).

    See Also
    --------
    lagged_ami

    References
    ----------
    Arbanel, H. D. (1996). *Analysis of Observed Chaotic Data* (p. 28). New York: Springer.

    """
    x, y = _vector_pair(x, y)
    if x.shape[0] != y.shape[0]:
        raise ValueError('timeseries must have the same length')

    return metrics.mutual_info_score(None, None, contingency=np.histogram2d(x, y, bins=n_bins)[0])


def lagged_ami(x, min_lag=0, max_lag=None, lag_step=1, n_bins=10):
    """Calculate the average mutual information between $x(t)$ and $x(t + \tau)$, at multiple values of $\tau$.

    Parameters
    ----------
    x : array-like
        $x(t)$.
    min_lag : int, optional
        The shortest lag to evaluate, in units of the sampling period $h$ of $x(t)$.
    max_lag : int, optional
        The longest lag to evaluate, in units of $h$.
    lag_step : int, optional
        The step between lags to evaluate, in units of $h$.
    n_bins : int
        The number of bins to use when computing the joint histogram in order to calculate mutual information.
        See |ami|.

    Returns
    -------
    lags : ndarray
        The evaluated lags $\tau_i$, in units of $h$.
    amis : ndarray
        The average mutual information between $x(t)$ and $x(t + \tau_i)$.

    See Also
    --------
    ami

    """
    if max_lag is None:
        max_lag = x.shape[0]//2
    lags = np.arange(min_lag, max_lag, lag_step)

    amis = [ami(reconstruct(x, lag, 2), n_bins=n_bins) for lag in lags]
    return lags, np.array(amis)


def _vector_pair(a, b):
    a = np.squeeze(a)
    if b is None:
        if a.ndim != 2 or a.shape[1] != 2:
            raise ValueError('with one input, array must have be 2D with two columns')
        a, b = a[:, 0], a[:, 1]
    return a, np.squeeze(b)


def _vector(x):
    x = np.squeeze(x)
    if x.ndim != 1:
        raise ValueError('x(t) must be a 1-dimensional signal')
    return x
