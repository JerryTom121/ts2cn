# code to perform stationarity analysis

# importing libraries
from numpy import array








def chi2_modified(x, y, dof=None):
	"""
##
# Perform the Chi Square test to compare two binned data
##
# x		Array. The frequencies, not density, for each bin. The vectors x and y should have the same support
# y		Array. The frequencies, not density, for each bin. The vectors x and y should have the same support
# dof	Number. Optional. Degree of freedom of the distribution
##
# Return a tuple with the statistic and its p-value
##
	"""
	from scipy.stats import chi2
	stat = sum((x - y)**2 / (x+y))
	if not dof:
		dof = len(x)
		if sum(x) == sum(y):
			dof = len(x) -1
	
	return (stat, 1-chi2.cdf(stat, df=dof))



def gp_integral_correlation(x, epsilon, thau_c=0, norm=None):
	"""
##
# Compute the Grassberger-Procaccia integral correlation
##
# x			Array. Phase space, each array element is an point in the phase space. The points need to be sorted regarding the trajectories.
# epsilon	Number. Distance epsilon
# thau_c	Number. An upper limit regarding the points to compute the integral correlation. 
#			Zero to compute the regular GP integral correlation. 
#			As suggested by Theiler, J. (1990). Estimating fractal dimension. JOSA A, 7(6), 1055-1073.,
#			see http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm for more details.
# norm		Function. A function that compute the norm between two points. Receives two points and return its distance
##
# Return a number, the integral correlation
# It is advised to have thau_c > thau*(2/N)^2/m, where thau is the 
[1] J. Theiler, “Spurious dimension from correlation algorithms applied to limited time-series data,” Phys. Rev. A, vol. 34, no. 3, pp. 2427–2432, 1986.
##

	"""
	from scipy.spatial.distance import pdist
	if not norm:
		norm = 'euclidean'
	# this mask is used to only compute the correlation for the given thau_c
	# it can be seem as only considering thau_c points above the diagonal
	N = len(x)
	mask = [[thau_c*[False], (N-i-thau_c)*[True]] for i in range(1, N-thau_c)]
	mask.append([i*[False] for i in range(thau_c, 0, -1)])
	# the if inside the listcomprehesion is needed to avoid try to hstack empty lists
	mask = hstack([hstack(i) for i in mask if i]) == 1
	c = sum(pdist(x, metric=norm)[mask] > epsilon) * 2/((N-thau_c)*(N-thau_c-1))
	
	return c


#  w=1; N=4;tt = [[w*[False], (N-i-w)*[True]] for i in range(1,N-w)]; tt.append([i*[False] for i in range(w, 0, -1)]); tt; yy=hstack([hstack(i) for i in tt]) == 1; yy; len(yy); squareform(yy)


def cross_correlation_integral(x, y, epsilon, norm=None):
	"""
##
# Compute the cross-correlation integral, it is a modified form of the Grassberger-Procaccia integral correlation
##
# x			Array. Phase space, each array element is an point in the phase space. The points need to be sorted regarding the trajectories.
# y			Array. Phase space, each array element is an point in the phase space. The points need to be sorted regarding the trajectories.
# epsilon	Number. Distance epsilon
# norm		Function. A function that compute the norm between two points. Receives two points and return its distance
##
# Return a number, the integral correlation between two time series (windows)
# It is computed comparing every point of x with every point of y, points in its respective phase space
# In this case there is no parameter thau_c
##

	"""
	from scipy.spatial.distance import euclidean, cdist
	if not norm:
		norm = euclidean
	#c = sum([norm(i, j) < epsilon for i in x for j in y])
	N_x = len(x)
	N_y = len(y)
	
	c = sum(sum(cdist(x, y) < epsilon))
	
	return c/(N_x*N_y)



def segment_closeness(x, y, m, epsilon_a, epsilon_b, thau_a, thau_b, thau_c):
	"""
##
# Compute the closeness between two windows (time series segments). The dynamical measure of the closeness between two orbits
##
# x			Array. Time series segment
# y			Array. Time series segment
# m			Number. Largest embedding dimension to compute the closeness, but all dimension from 1 to m will be used
# epsilon	Array. List with the epsilons to be used in the cross correlation integral
# thau		Number. The time lag used in the phase space reconstruction
# thau_c	Number. NOT USED YET. An upper limit regarding the points to compute the integral correlation in the grassberger-procaccia algorithm. 
#			Zero to compute the regular GP integral correlation. 
#			As suggested by Theiler, J. (1990). Estimating fractal dimension. JOSA A, 7(6), 1055-1073.,
#			see http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm for more details.
##
# Return a 2-dimensional array the computed segment closeness (D) between the two windows
# The first dimension is related to the epsilon and the second to the embedding dimension (m) used
##

[1] R. Manuca and R. Savit, “Stationarity and nonstationarity in time series analysis,” Phys. D Nonlinear Phenom., vol. 99, no. 2–3, pp. 134–161, Dec. 1996.


imp.reload(st); t=date.now(); b=[st.segment_closeness(hstack(tss['O'][i]), hstack(tss['O'][j]), 10, (std(hstack(tss['O'][0:1]))*0.1)*array([1,1.5,2,2.5]), 1, 0) for i in range(3) for j in range(i+1, 3)]; b; print(date.now() -t);

	"""
	from ts2cn.ts.phase_space import reconstruct
	from numpy import sqrt, absolute, maximum, append, nan

	from datetime import datetime as date
	t = date.now()
	u = date.now()
	
	D = []
	# to improve the performance allowing to compute several epsilon here, if not for different epsilons will require
	# the same reconstruction
	for e in range(len(epsilon_a)):
		
		C_xx = [0]*(m+2)
		C_yy = [0]*(m+2)
		C_xy = [0]*(m+2)
		
		# determine the max embedding dimension, based on the amount of data available
		max_m = -1
		
		# iterate over m (the embedding dimension)
		for emb_dim in range(1, m+2):
			# Reconstructing the phase space for the given parameters
			phs_x = reconstruct(x, lag=thau_a, n_dims=emb_dim)
			phs_y = reconstruct(y, lag=thau_b, n_dims=emb_dim)
			
			# skip this there are not enough data to reconstruct
			if len(phs_x) == 0 or len(phs_y) == 0:
				C_xx[emb_dim] = nan
				C_yy[emb_dim] = nan
				C_xy[emb_dim] = nan
				if max_m != -1:
					max_m = emb_dim
				continue

			# Computing the correlation integral for both phase spaces
			# TODO maybe use the gp_integral_correlation for calculate C_xx and C_yy, to account for the thau_c
			#C_x = gp_integral_correlation(phs_x, epsilon, thau_c, norm=None)
			#C_y = gp_integral_correlation(phs_y, epsilon, thau_c, norm=None)
			C_xx[emb_dim] = cross_correlation_integral(phs_x, phs_x, epsilon_a[e], norm=None)
			C_yy[emb_dim] = cross_correlation_integral(phs_y, phs_y, epsilon_b[e], norm=None)
			C_xy[emb_dim] = cross_correlation_integral(phs_x, phs_y, (epsilon_a[e]+epsilon_b[e])/2, norm=None)

		# As defined in the paper "Stationarity and nonstationarity in time series analysis, S_m = C_(m+1) / C_m
		# and R_m = S_m / S_(m-1), so R_m = C_(m+1) * C_m(m-1) / C_m^2
		# Because R_m = S_m / S_(m-1)
		# for m == 1 S_0 (predictability for m=0) is not defined so use C_x and C_y instead of R_1 (ratio of predictability for m=1)
		R_xx = array([0.0]*(m+1))
		R_xx[1] = C_xx[1]
		R_yy = array([0.0]*(m+1))
		R_yy[1] = C_yy[1]
		R_xy = array([0.0]*(m+1))
		R_xy[1] = C_xy[1]
		
		# TODO is it possible to transform these operations into a vector form (numpy array)?
		# I think it's not necessary, usually m is not a large number so the dimension of R is short
		for emb_dim in range(2, m+1):
			R_xx[emb_dim] = C_xx[emb_dim-1] * C_xx[emb_dim+1] / (C_xx[emb_dim]*C_xx[emb_dim])
			R_yy[emb_dim] = C_yy[emb_dim-1] * C_yy[emb_dim+1] / (C_yy[emb_dim]*C_yy[emb_dim])
			R_xy[emb_dim] = C_xy[emb_dim-1] * C_xy[emb_dim+1] / (C_xy[emb_dim]*C_xy[emb_dim])
		
		# Now compute the closeness between segments based on the ratio of predictability R_m
		D_xy = sqrt(absolute( R_xx - R_xy ))
		D_yx = sqrt(absolute( R_yy - R_xy ))
		
		# check if the maximum embedding dimension was reached
		# then set D_xy[max_m:] = -1 and D_yx[max_m:] = -1
		if max_m != -1:
			D_xy[max_m:] = -1
			D_yx[max_m:] = -1
		
		D.append(maximum(D_xy, D_yx))
	return array(D)

"""
imp.reload(st); t=date.now(); st.segment_closeness(hstack(tss['O'][0:4]), hstack(tss['O'][4:8]), 10, std(hstack(tss['O'][0:8]))*0.1, 1, 0); print(date.now() -t)
"""



def all_segment_closeness(ts, w_size, m, epsilon, thau, thau_c):
	"""
##
# Compute the segment closeness between all window pairs
##
# ts		Array. The Time series
# w_size	Number. The window size of the time series segments
# m			Number. Largest embedding dimension to compute the closeness, but all dimension from 1 to m will be used
# epsilon	Array. List with the epsilons to be used in the cross correlation integral
# thau		Array. Number. The time lag used in the phase space reconstruction
# thau_c	Number. NOT USED YET. An upper limit regarding the points to compute the integral correlation in the grassberger-procaccia algorithm. 
#			Zero to compute the regular GP integral correlation. 
#			As suggested by Theiler, J. (1990). Estimating fractal dimension. JOSA A, 7(6), 1055-1073.,
#			see http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm for more details.
##
# Return a 2-dimensional array.
# The D matrix where each element ij represents the segment closeness betweeen windows i and j
# Each element is a 2-dimensional array in itself, d_k_l, with the first dimension representing the epsilon 
# and the second the embedding dimension used to reconstruct and compute the integral correlation. See
# function segment_closeness for more details
##
Reference
[1] R. Manuca and R. Savit, “Stationarity and nonstationarity in time series analysis,” Phys. D Nonlinear Phenom., vol. 99, no. 2–3, pp. 134–161, Dec. 1996.


imp.reload(st); t=date.now(); d=st.all_segment_closeness(hstack(tss['O'][0:3]), len(tss['O'][0]), 10, (std(hstack(tss['O'][0:3]))*0.1)*array([1,1.5,2,2.5]), 1, 0); d; print(date.now() -t);

	"""
	from ts2cn.ts.misc import ts_windowed_indices as w_indices
	from scipy.spatial.distance import squareform
	
	# get the window indices from the entire time series
	w_i = w_indices(len(ts), window_size=w_size, overlap=0)
	w_n = len(w_i)
	
	D = array([segment_closeness(array(ts[ w_i[i] : (w_i[i]+w_size) ]), array(ts[ w_i[j] : (w_i[j]+w_size) ]), m, epsilon[:, i], epsilon[:, j], thau[i], thau[j], thau_c) 
				for i in range(w_n) for j in range(i+1, w_n)])
	return array([[squareform(D[:, e, emb_dim]) for emb_dim in range(m+1)] for e in range(len(epsilon))])




def dynamic_similarity(D, ths=1.414):
	"""
##
# Compute the dynamic similatiry between all pair of windows
##
# D		4-dimensional array as returned by the function all_segment_closeness. See function all_segment_closeness for more details
# ths	Number. The threshold used to compare the dynamical similarity between two windows. Default it 1.414
##
# Return a 2-dimensional array.
# The M matrix, dynamical similarity, where each element ij represents the dynamical similarity between windows i and j
# Each element is a 2-dimensional array in itself, m_k_l, with the first dimension representing the epsilon
# and the second the embedding dimension used to reconstruct and compute the segment closeness. See function
# all_segments_closeness for more details
##

Reference
[1] R. Manuca and R. Savit, “Stationarity and nonstationarity in time series analysis,” Phys. D Nonlinear Phenom., vol. 99, no. 2–3, pp. 134–161, Dec. 1996.

imp.reload(st); t=date.now(); e=st.dynamic_similarity(d); e; print(date.now() -t);
	"""
	from ts2cn.ts.phase_space import reconstruct
	from numpy import sqrt, ndarray, array, any, nan, isnan
	from scipy.spatial.distance import squareform
	
	# see reference for explaination of the value sqrt(2)
	# instead of using the fixed value sqrt(2) now it is a parameter, ths, defined by the user
	# empirically I noticed that for my application domain the usage of ths=1.414 is too permissive
	
	M = []
	
	# check for not a number, nan, in D and return a matrix filled with -1 as a marker
	nan_matrix = ndarray(shape=(len(D[0][0]), len(D[0][0])), dtype=float, buffer=array([-1.]*len(D[0][0])*len(D[0][0])))
	
	M = array([
		[ squareform([ sum( (D[e, m, :, i] - D[e, m, :, j]) > ths) for i in range(len(D[e, m])) for j in range(i+1, len(D[e, m])) ])
		 
			if not any(isnan(D[e, m])) else nan_matrix for m in range(len(D[e]))] 
				for e in range(len(D))
		])
	return M



def recurrence_analysis(ts, w_size, m, epsilon, thau, thau_c, ths=1.414):
	"""
##
# Perform the recurrence analysis on the time series
##
# ts		Array. The Time series
# w_size	Number. The window size of the time series segments
# m			Number. Largest embedding dimension to compute the closeness, but all dimension from 1 to m will be used
# epsilon	Array. List with the epsilons to be used in the cross correlation integral
# thau		Number. The time lag used in the phase space reconstruction
# thau_c	Number. NOT USED YET. An upper limit regarding the points to compute the integral correlation in the grassberger-procaccia algorithm. 
#			Zero to compute the regular GP integral correlation. 
#			As suggested by Theiler, J. (1990). Estimating fractal dimension. JOSA A, 7(6), 1055-1073.,
#			see http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm for more details.
# ths		Number. The threshold used to compare the dynamical similarity between two windows. Default it 1.414
##
# Return a 3 element tuple. The first element is the matrix D, the second the matrix M, and the third the clustering results for every setting.
# Each matrix and the clustering follows the index order as the D and M matrix, first dimension regarding thau, second w_size, third ths, fourth
# epsilon and fifth the embedding dimension m.
#
# For details of D and M matrix see all_segments_closeness and dynamic_similarity functions respectively.
#
# For a given setting, thau, w_size, ths, epsilon and m, the clustering structure gives the cluster index corresponding to the respective window
# for instance, cl[0][0][0][0][0] = [0, 0, 1, 1] means that the first and second windows are in the cluster zero and the third and fourth in the 
# cluster one.
# 
##

Reference
[1] R. Manuca and R. Savit, “Stationarity and nonstationarity in time series analysis,” Phys. D Nonlinear Phenom., vol. 99, no. 2–3, pp. 134–161, Dec. 1996.

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(tss['O'][0:3]), 2048, 10, (std(hstack(tss['O'][0:3]))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:03:09.604900

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(tss['O'][0:3]), 1024, 10, (std(hstack(tss['O'][0:3]))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:03:47.005395

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(tss['O'][0:3]), 512, 10, (std(hstack(tss['O'][0:3]))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:04:30.933287

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(tss['O'][0:3]), 256, 10, (std(hstack(tss['O'][0:3]))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:05:28.613648

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(tss['O'][0:3]), 128, 10, (std(hstack(tss['O'][0:3]))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:09:23.423130

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(tss['O'][0:3]), 64, 10, (std(hstack(tss['O'][0:3]))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:26:58.769833



imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))), 2048, 10, (std(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))))*0.1)*array([0.7, 0.1]), 1, 0); r; print(date.now() -t);
0:07:02.070001

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))), 1024, 10, (std(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))))*0.1)*array([0.7, 0.1]), 1, 0); r; print(date.now() -t);
0:07:56.598315

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))), 512, 10, (std(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))))*0.1)*array([0.7, 0.1]), 1, 0); r; print(date.now() -t);
0:09:21.180265


imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))), 2048, 10, (std(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:14:59.402063

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))), 1024, 10, (std(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);
0:16:24.618490

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))), 512, 10, (std(hstack(vstack((tss['O'][0:3], tss['S'][0:3]))))*0.1)*array([0.7, 0.5, 0.2, 0.1]), 1, 0); r; print(date.now() -t);





imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack((ts1[0], ts2[0])), 1024, 10, (std(hstack((ts1[0], ts2[0])))*0.1)*array([20, 10, 5, 2]), 1, 0); 'r'; print(date.now() -t); sum(sum(r[0])); sum(sum(r[1])); sum(sum(r[2]));

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack((ts1[0], ts2[0])), 128, 10, (std(hstack((ts1[0], ts2[0])))*0.1)*array([20, 10, 5, 2]), 1, 0); 'r'; print(date.now() -t); sum(sum(r[0])); wrong=[any(isnan(i)) for i in r[0]]; wrong; sum(sum(r[1])); sum(sum(r[2])); len(sum(sum(r[2]))); sum(sum(r[2][where(wrong == False)]));



imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack((ts1[0], ts2[0])), 128, 10, (std(hstack((ts1[0], ts2[0])))*0.1)*array([20, 10, 5, 2]), 1, 0); 'r'; print(date.now() -t); sum(sum(r[0])); wrong=array([any(isnan(i)) for i in r[0]]); wrong; sum(sum(r[1])); sum(sum(r[2])); len(sum(sum(r[2]))); sum(r[2][logical_not(wrong)])

imp.reload(st); t=date.now(); r=st.recurrence_analysis(hstack((ts1[0], ts2[0])), 128, 10, (std(hstack((ts1[0], ts2[0])))*0.1)*array([20, 10, 5, 2]), 1, 0); 'r'; print(date.now() -t); sum(sum(r[0])); wrong=array([any(isnan(i)) for i in r[0]]); wrong; sum(sum(r[1])); sum(sum(r[2])); len(sum(sum(r[2]))); sum(sum(r[2][logical_not(wrong)]));
	"""
	from ts2cn.ts.phase_space import reconstruct
	from numpy import sqrt, average, where, transpose
	
	# compute the matrix of closeness between windows
	D = all_segment_closeness(ts, w_size, m, epsilon, thau, thau_c)
	
	# compute the matrix of dynamic similarity between windows
	M = dynamic_similarity(D, ths)
	
	clusters = []
	# iterate over m and epsilon to navigate through M
	for e in range(len(M)):
		clusters.append([])
		
		for emb_dim in range(len(M[0])):
		
			# I tried to follow the same nomenclature as the article says, see pg 8 of the reference

			# compute the average of M for the current epsilon and embedding dimension
			M_avg = sum(M[e, emb_dim]) / len(M[e, emb_dim])
			
			# initializing variable of the recurrence analysis
			G_m = []
			n_windows = len(M[e, emb_dim])
			clusters[e].append( [-1]*n_windows )
			sqrt_n = sqrt(n_windows)
			# m represent the current cluster being analyzed, and it starts at 0
			cl = -1
			# creating list of unclustered windows
			n_u = list(range(n_windows))
			
			while len(n_u):
				
				# step 2 in the recurrence analysis, see reference pg 8
				# compute the close neighbors, N_i, for every unclustered window
				N_i = [sum(sqrt_n > M[e, emb_dim, i, n_u]) for i in n_u]
				
				# step 3 in the recurrence analysis
				# increase the number of clusters, get the window with the greatest N_i
				# and put it into the current cluster, cl
				G_m.append([])
				cl += 1
				i = N_i.index(max(N_i))
				G_m[cl].append( n_u[i] )
				n_u.pop(i)
				
				# checking if there are still unclustered windows
				if not n_u:
					break
				
				# step 4 in the recurrence analysis
				# compute the average of M for each window
				# get the windows i consistent with the relation
				# <M_j_i> W_j in G_m[cl] < sqrt(n)
				# and put them into the current cluster m
				# I think this process has to be iterative, i.e, compute
				# one window at a time because this changes the average <M_j_i>
				# with the addition of a window into the clustered windows array
				
				# TODO check again and garantee that the 4-th step is correct, I still have doubts
				M_avg = sum(transpose(M[e, emb_dim, n_u, :][:, G_m[cl]])) / len(G_m[cl])
				min_M_avg = min(M_avg)
				i_min_M_avg = where(M_avg == min_M_avg)[0][0]
				
				while (sum(M[e, emb_dim, G_m[cl], i_min_M_avg]) / len(G_m[cl])) < sqrt_n:
					# if the window with the index i_min_M_avg respects the condition
					# then insert it into the m-th cluster and try again with the
					# unclustered windows
					G_m[cl].append( n_u[i_min_M_avg] )
					n_u.pop(i_min_M_avg)
					
					# checking if there are still unclustered windows
					if not n_u:
						break
					M_avg = sum(transpose(M[e, emb_dim, n_u, :][:, G_m[cl]])) / len(G_m[cl])
					min_M_avg = min(M_avg)
					i_min_M_avg = where(M_avg == min_M_avg)[0][0]
				
				# TODO the doubt is regarding how to compute M_avg
				"""
				older version of the step 4
				min_M_avg = min(M_avg[n_u])
				i_min_M_avg = where(M_avg == min_M_avg)[0][0]
				while (sum(M[e, emb_dim, clusters[e][emb_dim][cl], i_min_M_avg])/len(clusters[e][emb_dim][cl])) < sqrt_n:
					# if the window with the index i_min_M_avg respects the condition
					# then insert it into the m-th cluster and try again with the
					# unclustered windows
					n_u.pop(i_min_M_avg)
					clusters[e][emb_dim][cl].append(i_min_M_avg)
					min_M_avg = min(M_avg[n_u])
					i_min_M_avg = where(M_avg == min_M_avg)[0][0]
				"""
				
				# step 5 is the repetition of steps 2-4 until all windows are clustered
				# hence the first while loop
			
			# setting the clusters index to the windows positions
			for i in range(len(G_m)):
				for j in G_m[i]:
					clusters[e][emb_dim][j] = i

	return (D, M, array(clusters))





def select_recur_parameters(tss, grd_truth, w_sizes, m, epsilon, thau={}, thau_c=0, ths=[1.414, 1, 0.5, 0.25, 0.1, 0.01]):
	"""
##
# Evaluate the clustering quality of the recurrence analysis in order to determine a set of parameters for it
##
# TODO fix documentation
# tss		A 2-dimensional array. Each element, line, of the array is an array with two time series glued together 
#			each pertaining to a different cluster, they are divided in the middle.
# w_sizes	Array. An array with the window sizes to analyse
# m			Number. Largest embedding dimension to compute the closeness, but all dimension from 1 to m will be used
# epsilon	A 2-dimensionla array. Each element of the array, line, is a list with the epsilons to be used in the 
#			cross correlation integral for the respective element in the tss array
# thau		A 2-dimensional array. Each element, line, is a list of time lags used in the phase space reconstruction
#			for the respective element in the tss array
# thau_c	Number. NOT USED YET. An upper limit regarding the points to compute the integral correlation in the grassberger-procaccia algorithm. 
#			Zero to compute the regular GP integral correlation. 
#			As suggested by Theiler, J. (1990). Estimating fractal dimension. JOSA A, 7(6), 1055-1073.,
#			see http://www.scholarpedia.org/article/Grassberger-Procaccia_algorithm for more details.
# ths		Array. List of thresholds to analyse, see dynamic_similarity function for more details
# grd_truth			List. A list the to ground truth labels for each instance of the time series
#
# - see recurrence_analysis function for more details about the parameters
##
# Return 
##

Reference
[1] R. Manuca and R. Savit, “Stationarity and nonstationarity in time series analysis,” Phys. D Nonlinear Phenom., vol. 99, no. 2–3, pp. 134–161, Dec. 1996.

Usage example:

from numpy import array, std, average, hstack, logical_not, isnan, nan, any, where
from datetime import datetime as date
import imp
from datasets import bonn
from ts2cn.ts import stationarity as st
from ts2cn.ts.phase_space import compute_lag

tss = bonn.load('/home/daniel/docs/mestrado/dataset/epilepsy/quali/bonn/data/')
tss_orig = tss.copy()
#tss = [hstack((tss['O'][i], tss['S'][i])) for i in range(30)]
# server
tss = [hstack((tss['O'][i], tss['S'][i])) for i in range(30)]
w_sizes = [128, 256, 512, 1024, 2048, 4096]
ths = [1.414, 1, 0.5, 0.25, 0.1, 0.01]
m = 20
epsilon = [std(ts)*0.1*array([20, 10, 5, 2, 1]) for ts in tss]
thau = [[1, compute_lag(ts)] for ts in tss]

first_ts = 'O'
second_ts = 'S'
label_tmp = hstack( ([first_ts]*len(tss_orig[first_ts][0]), [second_ts]*len(tss_orig[second_ts][0])) )
labels = [label_tmp for i in tss]

imp.reload(st); t=date.now(); r=st.select_recur_parameters(tss[:1], labels[:1], w_sizes[-3:], 10, epsilon[-4:], thau[:2], ths=ths[:1], thau_c=0); print(date.now() -t);

# server
imp.reload(st); t=date.now(); r_128_1=st.select_recur_parameters(tss, labels, w_sizes[0:1], m, epsilon, thau, ths=ths, thau_c=0); print(date.now() -t);
imp.reload(st); t=date.now(); r_128_2=st.select_recur_parameters(tss, labels, w_sizes[0:1], m, epsilon, thau, ths=ths, thau_c=0); print(date.now() -t);
imp.reload(st); t=date.now(); r_128_3=st.select_recur_parameters(tss, labels, w_sizes[0:1], m, epsilon, thau, ths=ths, thau_c=0); print(date.now() -t);
imp.reload(st); t=date.now(); r_256=st.select_recur_parameters(tss, labels, w_sizes[1:2], m, epsilon, thau, ths=ths, thau_c=0); print(date.now() -t);
imp.reload(st); t=date.now(); r_512=st.select_recur_parameters(tss, labels, w_sizes[2:3], m, epsilon, thau, ths=ths, thau_c=0); print(date.now() -t);
imp.reload(st); t=date.now(); r_1024=st.select_recur_parameters(tss, labels, w_sizes[3:], m, epsilon, thau, ths=ths, thau_c=0); print(date.now() -t);
	"""
	from sklearn.metrics import adjusted_rand_score as rand_id
	from numpy import std, average, array, hstack, isnan, nan, any, where, ndarray
	from math import ceil

	from ts2cn.ts import misc
	import imp
	imp.reload(misc)
	
	rs_index = []
	rs_index_perf = []
	rs_index_perf_sd = []
	rs_index_avg = []
	rs_index_sd = []
	rs_clusters = []
	rs_clusters_avg = []
	rs_D = []
	rs_M = []
	for j in range(len(thau[list(thau.keys())[0]][0][0])):
		
		rs_index.append([])
		rs_index_perf.append([])
		rs_index_perf_sd.append([])
		rs_index_avg.append([])
		rs_index_sd.append([])
		rs_clusters.append([])
		rs_clusters_avg.append([])
		rs_D.append([])
		rs_M.append([])
		
		for w_size_i, w_size in enumerate(w_sizes):
			
			rs_index[j].append([])
			rs_index_perf[j].append([])
			rs_index_perf_sd[j].append([])
			rs_index_avg[j].append([])
			rs_index_sd[j].append([])
			rs_clusters[j].append([])
			rs_clusters_avg[j].append([])
			rs_D[j].append([])
			rs_M[j].append([])
			
			for dyn_ths_i, dyn_ths in enumerate(ths):
			
				rs_index[j][w_size_i].append( [] )
				rs_index_perf[j][w_size_i].append( [] )
				rs_index_perf_sd[j][w_size_i].append( [] )
				rs_index_avg[j][w_size_i].append( array([[0.]*(m+1)]*len(epsilon[0])) )
				rs_index_sd[j][w_size_i].append( [] )
				rs_clusters[j][w_size_i].append( [] )
				rs_clusters_avg[j][w_size_i].append( array([[[0.]*int(len(tss[0])/w_size)]*(m+1)]*len(epsilon[0])) )
				
				# create a zero matrix for averaging D and M, closeness and dynamical similarity matrices
				number_w = int(len(tss[0])/w_size)
				zero_matrix = ndarray(shape=(number_w, number_w), dtype=float, buffer=array([0.]*number_w*number_w))
				
				rs_D[j][w_size_i].append( array([[ zero_matrix ]*(m+1)]*len(epsilon[0])) )
				rs_M[j][w_size_i].append( array([[ zero_matrix ]*(m+1)]*len(epsilon[0])) )
				
				cl_index_sd = []
				cl_index_perf_sd = []
				
				for i, ts in enumerate(tss):
					
					print("(ts, len(ts), w_size, m, i, epsilon[w_size_i][:, i], (i,j), array(thau[w_size])[i, :, j], thau_c, dyn_ths)")
#					print(epsilon[w_size][i])
#					print(array(thau[w_size])[i, :, j])
					print((ts, len(ts), w_size, m, i, epsilon[w_size_i][:, i], (i,j), array(thau[w_size])[i, :, j], thau_c, dyn_ths))
					
#					D, M, cl = recurrence_analysis(ts, w_size, m, epsilon[i], array(thau[w_size])[i, :, j], thau_c, ths=dyn_ths)
					print(epsilon[w_size_i].shape)
					D, M, cl = recurrence_analysis(ts, w_size, m, epsilon[w_size_i][:, i], array(thau[w_size])[i, :, j], thau_c, ths=dyn_ths)
					size_cl = ceil(len(cl[0][0])/2)

					# grd_truth passed as parameter
					# given the ground truth for each element of the time series, compute the mode, most common element, for the windows
					grd_truth_ts = misc.w_labels(grd_truth[i], w_size)
#					grd_truth_ts = [ max( set(grd_truth[i][ (k*w_size):((k+1)*w_size) ]), key=list(grd_truth[i][ (k*w_size):((k+1)*w_size) ]).count )
#										for k in range(int(len(ts)/w_size))]

					# evaluate the clustering performance with the adjusted rand index
					cl_index = array([[rand_id(grd_truth_ts, cl_e_m) for cl_e_m in cl_e] for cl_e in cl])
					cl_index_perf = array([[ cluster_eval(cl_e_m, grd_truth_ts) for cl_e_m in cl_e] for cl_e in cl])
					print(('cl_index', cl_index))
					# this sum operation is allowed because these elements are array and element-wise
					# summation are defined
					rs_index[j][w_size_i][dyn_ths_i].append( cl_index )
					rs_index_avg[j][w_size_i][dyn_ths_i] += cl_index
					cl_index_sd.append( cl_index.copy() )
					cl_index_perf_sd.append( cl_index_perf.copy() )
#					rs_index_sd[j][w_size_i][dyn_ths_i].append( cl_index )
					rs_clusters[j][w_size_i][dyn_ths_i].append( cl )
					rs_clusters_avg[j][w_size_i][dyn_ths_i] += cl
					rs_D[j][w_size_i][dyn_ths_i] += D
					rs_M[j][w_size_i][dyn_ths_i] += M

					# TODO check for the necessity of this check, the cluster index in this case will be zero
					# find if there is any not a number, nan, in D
					#nan_indices = where(isnan(D))
					#nan_indices = set([(nan_indices[0][i], nan_indices[1][i]) for i in range(len(nan_indices[0]))])
				
				# take the average over all clusters
				rs_index_avg[j][w_size_i][dyn_ths_i] /= len(tss)
				print(('cl_index_sd FIM', cl_index_sd))
				cl_index_sd = array(cl_index_sd)
				rs_index_sd[j][w_size_i][dyn_ths_i] = [[std(cl_index_sd[:, cl_e, cl_e_m]) for cl_e_m in range(len(cl_index_sd[0][0]))] for cl_e in range(len(cl_index_sd[0])) ]
				cl_index_perf_sd = array(cl_index_perf_sd)
				rs_index_perf[j][w_size_i][dyn_ths_i] = [[average(cl_index_perf_sd[:, cl_e, cl_e_m]) for cl_e_m in range(len(cl_index_perf_sd[0][0]))] for cl_e in range(len(cl_index_perf_sd[0])) ]
				rs_index_perf_sd[j][w_size_i][dyn_ths_i] = [[std(cl_index_perf_sd[:, cl_e, cl_e_m]) for cl_e_m in range(len(cl_index_perf_sd[0][0]))] for cl_e in range(len(cl_index_perf_sd[0])) ]
				print(('rs_index_sd[][][] ', rs_index_sd[j][w_size_i][dyn_ths_i]))
				rs_clusters_avg[j][w_size_i][dyn_ths_i] /= len(tss)
				rs_D[j][w_size_i][dyn_ths_i] /= len(tss)
				rs_M[j][w_size_i][dyn_ths_i] /= len(tss)
				
				# turn them back to list
				rs_clusters_avg[j][w_size_i][dyn_ths_i] = list([list(k) for k in rs_clusters_avg[j][w_size_i][dyn_ths_i]])
				rs_D[j][w_size_i][dyn_ths_i] = list([list(k) for k in rs_D[j][w_size_i][dyn_ths_i]])
				rs_M[j][w_size_i][dyn_ths_i] = list([list(k) for k in rs_M[j][w_size_i][dyn_ths_i]])
					
	rs_index_avg = array(rs_index_avg)
	rs_index_sd = array(rs_index_sd)
	rs_clusters_avg = array(rs_clusters_avg)
	rs_D = array(rs_D)
	rs_M = array(rs_M)
	rs_index_perf = array(rs_index_perf)
	rs_index_perf_sd = array(rs_index_perf_sd)
	

	# ATTENTION the usage of dtype=object allows the incorporation of several lists with different shapes
	# into a single array but this array loses its ability of understading its elements
	# like, the use of where and isnan for this array is no longer possible, it's needed to iterate over this
	# array in order to use this methods
	# I used this approach so I can use a single save method, from numpy, to save this array into a file
	return (	rs_index_avg,
				rs_clusters_avg,				# meanless
				rs_D,
				rs_M,
				rs_index_sd,
				rs_clusters,
				rs_index,
				rs_index_perf,
				rs_index_perf_sd,
			)
	return array((	rs_index.tolist(),
					rs_clusters.tolist(),
					rs_D.tolist(),
					rs_M.tolist(),
					rs_index_sd.tolist(),
			), dtype=object)

	 
	#array([[cl_index(j, grd2) for j in i] for i in r[2]])[1,20]
	#r=st.recurrence_analysis(hstack((ts1[0], ts2[0])), 2048, 20, (std(hstack((ts1[0], ts2[0])))*0.1)*array([20, 10, 5, 2]), 1, 0, ths=0.01); 
	#'r'; print(date.now() -t); sum(sum(r[0])); wrong=array([any(isnan(i)) for i in r[0]]); wrong; sum(sum(r[1])); sum(sum(r[2])); len(sum(sum(r[2]))); sum(sum(r[2][logical_not(wrong)]))



def cluster_eval(cluster, labels):
	"""
##
# Return the accuracy of the clusters with the given labels
##
# cluster:	List.
# labels:	List.
##
# Returns a float between 0 and 1 as a accuracy measure.
# Since the clustering process can return more than the number os labels, this function substitutes
# the cluster value for the label that appears the most for this given cluster
# example:
# for this cluster result [0, 0, 1, 2, 3, 3, 2, 2] and corresponding label ['no', 'no', 'sz', 'sz', 'no', 'no', 'sz', 'sz']
# the function translate the clustering results to ['no', 'no', 'sz', 'sz', 'no', 'no', 'sz', 'sz'] and then compares with the labels
#
# in another example, the clustering [1, 1, 0, 0, 1, 1, 0, 1] for the same label ['no', 'no', 'sz', 'sz', 'no', 'no', 'sz', 'sz'],
# is translated to ['no', 'no', 'sz', 'sz', 'no', 'no', 'sz', 'no'], here only the last element is wrongly clustered
# 
	"""
	cl_dict = dict([(cl, []) for cl in cluster])
	for i, cl in enumerate(cluster):
		cl_dict[cl].append( labels[i] )
	
	for cl in cl_dict:
		cl_dict[cl] = max( set(cl_dict[cl]), key=list(cl_dict[cl]).count )
	
	cl_labels = [cl_dict[cl] for cl in cluster]
	return sum(array(cl_labels) == array(labels)) / len(cluster)





def generic_stat_test(ts_data, alphas, labels, target_label='sz', id_data='N', save_dir='/tmp/st/', n_bins=100, w_sizes=[], 
						plot_title='', plot=False, only_chi2=True, target_str='Seizure only'):
	"""
##
# Perform a statistical test for stationarity
##
# ts_data:		
# id_data:		
# save_dir:		
# n_bins:		
# w_sizes:		
# plot:			
##
# Perform the stationary test based on the paper Stationarity of the EEG series - Blanco et at. 1995, but with a chi squared test and other modifications
##





	# TODO adjust the plot dimensions, width and height
	# print the time series and save in a file
	png(paste(save_dir, 'ts_', sub('/', '', id_data), '_', i, '.png', sep=''))
	plot(ts_i, type='l', main=paste('Time Series Type', sub('/', '', id_data), '- Number', i))
	dev.off()
	
	# TODO change this approach
	#	first compute all histograms, and possibly save then in a file, and then save the pictures of it
	
	# TODO perform some kind of statistical test with the histograms
	#	two consecutives histograms at a time to get the "zones" described in the paper
	
	# put some label and save in a file
	for(w_size in seq(from=round(min_window*ts_size), to=round(max_window*ts_size), by=round(step_window*ts_size)) ){
		for(j in 1:(round(ts_size/w_size)-1)){
			png(paste(save_dir, 'ts_', sub('/', '', id_data), '_', i, '_', w_size, '_', j, '.png', sep=''))
			hist(ts_i[(j*w_size):((j+1)*w_size)], breaks=n_bins)
			dev.off()
		}
	}
	"""

	from numpy import histogram
	from numpy import average, std, array, vstack, linspace, isnan
	from math import sqrt
	from scipy.stats import probplot, shapiro, anderson, kstest, normaltest, ks_2samp, entropy, chisquare, chi2
	from matplotlib import pyplot as plt
	from os import system

	from ts2cn.ts import misc
	
	# in case the time series passed are emnpty
	if not len(ts_data) or not len(ts_data[0]):
		return { 'w_sizes':[], 'hists': [], 'performance':{'ts':[], 'total':[], 'seizure':[]} }

	ts_size = len(ts_data[0])
	results = {
		'w_sizes': w_sizes,
		'hists': [],
		'performance': {'ts':[], 'total':[], 'seizure':[]}
	}

	import imp
	imp.reload(misc)
	imp.reload(plt)

	
	
	# iterate over ts_data with the var name ts_i
	for i in range(len(ts_data)):
		ts_i = ts_data[i]
		ts_size = len(ts_i)
		
		if plot:
			plt.figure(figsize=(27, 8))
			plt.plot(ts_i)
			title = plot_title + "%s - ts: %d" % (id_data, i)
			plt.title(title)
			if save_dir:
				plt.savefig("%sts_%s_%d.png" % (save_dir, id_data, i))
				plt.close()
			else:
				plt.show()

		# compute the histogram
		hists = []
		for w_i, w_size in enumerate(w_sizes):
			hists.append([])
			w_index = len(hists) -1
			power_w = []
			for j in range(0, int(ts_size/w_size)):
				ts = ts_i[ (j*w_size):((j+1)*w_size) ]
				if plot:
					plt.figure(figsize=(8, 19))
					plt.subplot(211)
					plt.hist(ts, bins=n_bins[w_i], normed=True)
				hist_w_size = histogram(ts, bins=n_bins[w_i], normed=True)
				w_avg = average(ts)
				w_std = std(ts)

				# compute normality tests
				test_chi2 = normaltest(ts)
				test_shapiro = shapiro(ts)
				stat_and = anderson(ts, dist='norm')
				test_anderson = (stat_and[0], stat_and[2][sum(stat_and[0] > stat_and[1]) -1]/100)
				test_ks = kstest(ts, cdf='norm')
				
				# measuring differences between the current and previous histograms
				kl_div = -1
				ks_2samp_test = (-1, -1)
				Dl1 = -1
				Dl2 = -1
				intersection = -1
				x2_1 = (-1, -1)
				x2_2 = (-1, -1)
				x2 = (-1, -1)
				eff_size = -1
				power = [-1] * len(alphas)
				if j > 0:
					# performing a Kolmogorov-Smirnov two samples test
					ts_prev = ts_i[ ((j-1)*w_size):(j*w_size) ]
					min_bin = min(min(ts_prev), min(ts))
					max_bin = max(max(ts_prev), max(ts))
					hist_bins = linspace(min_bin, max_bin, n_bins[w_i])
					
					# to only perform the analysis based on the Chi-Square statistic
					if not only_chi2:
						ks_2samp_test = ks_2samp(ts_prev, ts)
						
						# Calculating the Kullback-Leibler divergence
						# first it's needed to put histograms in a common bin scale, the same support
						hist_curr = histogram(ts, bins=hist_bins, normed=True)
						hist_prev = histogram(ts_prev, bins=hist_bins, normed=True)
						kl_div = entropy(hist_prev[0], hist_curr[0])
						
						# computing other distances, D_l1, D_l2, intersection and X^2
						diff_hists = [hist_prev[0][k] - hist_curr[0][k] for k in range(len(hist_prev[0]))]
						bin_width = hist_bins[1] - hist_bins[0]
						Dl1 = sum([abs(k) for k in diff_hists])*bin_width
						Dl2 = sqrt(sum([k**2 for k in diff_hists]))*bin_width
						intersection = sum([min([hist_prev[0][k], hist_curr[0][k]]) 
									for k in range(len(hist_prev[0]))])*bin_width
						x2_1 = chisquare([hist_curr[0][k]*bin_width for k in range(len(hist_curr[0])) if hist_prev[0][k] != 0],
								[hist_prev[0][k]*bin_width for k in range(len(hist_prev[0])) if hist_prev[0][k] != 0])
						x2_2 = chisquare([hist_prev[0][k]*bin_width for k in range(len(hist_curr[0])) if hist_curr[0][k] != 0],
								[hist_curr[0][k]*bin_width for k in range(len(hist_curr[0])) if hist_curr[0][k] != 0])

					# compute the modified Chi-squared test to see the differences between the two windows
					x2 = chi2_modified( histogram(ts_prev, bins=hist_bins, normed=False)[0], histogram(ts, bins=hist_bins, normed=False)[0] )
					
					# compute the effect size
					n = len(ts)
					df = n_bins[w_i] -1
					eff_size = sqrt(x2[0]/(n*df)) if x2[0] > 0 else 0

					##
					# compute the power of the test, using the Cramer's V statistic as effect size measure
					##
					# chi_crit is the critical value for the given alpha and degree of freedom
					chi_crits = [chi2.ppf(q=alpha, df=df, loc=0, scale=1) for alpha in alphas]
					# noncentrality parameter, the shift caused by the effect size
					loc = chi2.cdf(x=(x2[0]/df), df=df)
					# compute the power for every statistical test
					power = [1-chi2.cdf(x=chi_crit, df=df, loc=loc) for chi_crit in chi_crits]
				
				# save current histogram
				hists[w_index].append( {
					# histogram measures
					'avg': w_avg, 'std': w_std, 'freq': hist_w_size[0], 'bins': hist_w_size[1],
					# normality tests
					'normality': {
						'test_chi2': test_chi2, 'test_shapiro': test_shapiro, 'test_anderson': test_anderson, 'test_ks': test_ks, 
					},
					'stationarity': {
					# histograms differences
						'kl_div': kl_div, 'ks_2samp': ks_2samp_test, 'Dl1': Dl1, 'Dl2': Dl2,	'intersection': intersection, 'x2_1': x2_1, 'x2_2': x2_2, 'x2': x2,
					},
					'alphas': alphas,
					'power': power,
					'eff_size': eff_size
				} )
				
				# TODO move these plotting procedures to a separated function

				if plot:
					plt.title("%s - ts: %d - window size: %d - window index: %d \n avg: %s - std: %s\n\nKolmogorov-Smirnov 2 sample %s p-value %s\nD_L1 %s D_L2 %s Intersection %s\nChi-Square expec_prev %s p-value %s\nChi-Square expec_curr %s p-value %s\nChi-2 Mod stat %s p-value %s" % 
								(id_data, i, w_size, j, round(w_avg,1), round(w_std,1), 
	#							round(kl_div),
								round(ks_2samp_test[0], 2), round(ks_2samp_test[1], 3),
								round(Dl1, 3), round(Dl2, 3), round(intersection, 3), 
								round(x2_1[0], 3), round(x2_1[1], 3),
								round(x2_2[0], 3), round(x2_2[1], 3),
								round(x2[0], 3), round(x2[1], 3) ))
					
					plt.subplot(212)
					probplot((ts-w_avg)/w_std, dist='norm', plot=plt)
					plt.subplots_adjust(hspace=0.25)
					plt.title("Chi^2 %s p-value %s\nShapiro-Wilk %s p-value %s\nAnderson-Darling %s p-value < %s\nKolmogorov-Smirnov %s p-value %s" % 
								(round(test_chi2[0], 2), round(test_chi2[1], 3),
								round(test_shapiro[0], 2), round(test_shapiro[1], 3),
								round(test_anderson[0], 2), round(test_anderson[1], 3), 
								round(test_ks[0], 2), round(test_ks[1], 3)))
					"""
					plt.title("Probability Plot\nShapiro-Wilk %s p-value %s\nAnderson-Darling %s p-value < %s\nKolmogorov-Smirnov %s p-value %s" % (
							round(test_shapiro[0], 1), round(test_shapiro[1], 1),
							round(test_anderson[0], 1), round(test_anderson[1], 1), 
							round(test_ks[0], 1), round(test_ks[1], 1)))
					"""
					if save_dir:
						plt.savefig("%sts_%s_%d_%d_%d.png" % (save_dir, id_data, i, w_size, j))
						plt.close()
					else:
						plt.show()
						plt.close()
			

			if plot:
				# plot average, std, effect size and power for the given alphas
				plt.subplots(figsize=(30, 6))
				plt.xlim((-0.5, len(hists[w_index])))
				
				
				plt.subplot(141)
				plt.xlim((-0.5, len(hists[w_index])))
				plt.plot([i['avg'] for i in hists[w_index]], 'o-')
				plt.xlabel('Window index\n')
				plt.ylabel("Mean Value")
				
				plt.subplot(142)
				plt.xlim((-0.5, len(hists[w_index])))
				plt.plot([i['std'] for i in hists[w_index]], 'o-')
				plt.xlabel('Window index')
				plt.ylabel("Standard Deviation")

				plt.subplot(143)
				plt.xlim((-0.5, len(hists[w_index])))
				plt.plot(range(1, len(hists[w_index])), [i['eff_size'] for i in hists[w_index]][1:], 'o-')
				plt.xlabel('Window index')
				plt.ylabel("Effect Size")
				plt.title("%s - ts: %d - window size: %d - alphas %s" % (id_data, i, w_size, '  '.join([str(alpha) for alpha in alphas])), y=1.05)

				plt.subplot(144)
				markers = ['o', '^', 'v', '>', '<', 's', '*', 'D', '+', '|']
				markers_size = [8, 7.8, 7.6, 7.4, 7.2, 7, 6.8, 6.6, 6.4, 6.2, 6]
				power_w = array([ j['power'] for j in hists[w_index] ])
				print(('power_w', power_w))

				plt.ylim((-0.5, 1.1))
				plt.xlim((-0.5, len(power_w)))
				for a_i, a in enumerate(alphas):
					plt.plot(range(1, len(hists[w_index])), power_w[1:, a_i], '-', marker=markers[a_i], markersize=markers_size[a_i])
#				plt.plot([j['power'] for i,j in enumerate(hists[w_index])], '-', marker=markers[:len(hists[w_index])])
				plt.legend(labels=['alpha %.3f'% i for i in alphas], fontsize='x-small', framealpha=0, loc='lower center')
				plt.xlabel('Window index')
				plt.ylabel("Power")

				plt.subplots_adjust(wspace=0.8, left=0.125, right=0.9)
				
				if save_dir:
					plt.savefig("%sts_%s_%d_%d_avg_std.png" % (save_dir, id_data, i, w_size))
					plt.close()
				else:
					plt.show()
			
			
				# TODO save and plot average and std for a given window size over all windows
				base_name = "%sts_%s_%d_%d" % (save_dir, id_data, i, w_size)
				filenames = ' '.join(["%s_%s.png"%(base_name, j) for j in range(0, int(ts_size/w_size))])
				filename_avg = "%sts_%s_%d_%d_avg_std.png" % (save_dir, id_data, i, w_size)
				filename_hists = "%sts_%s_%d_%d.png" % (save_dir, id_data, i, w_size)
				filename_ts = "%sts_%s_%d.png" % (save_dir, id_data, i)
				# glue together all histogram plot a given window size
				system("convert %s +append %sts_%s_%d_%d.png" % (filenames, save_dir, id_data, i, w_size))
				system("convert %s %s +append %s_" % (filename_ts, filename_avg, filename_avg))
				system("convert %s_ %s -append %sts_%s_%d_%d_all.png" % (filename_avg, filename_hists, save_dir, id_data, i, w_size))

		
		##
		# Compute the performance for this given time series, for all window sizes
		##
		ts_labels = array(labels[i])
		ts_perf = [ 
				[ 
					compute_stat_perf(pvalues=[h['stationarity']['x2'][1] for h in hists[w_index]], target_label=target_label, grd_truth=ts_labels, w_size=w_size, alpha=alpha)
							for alpha_i, alpha in enumerate(alphas) 
				] 
								for w_index, w_size in enumerate(w_sizes) 
			]
		
		results['performance']['ts'].append( ts_perf )
		print();print()
		print(('w_sizes', w_sizes))
		print(('alphas', alphas))
		print(('ts_perf', ts_perf)); print()
		print(("results['performance']['ts']", results['performance']['ts'])); print(); print()
		# plot the performance measures for the all parameters on this particular time series
		if plot:
			for w_i, w_size in enumerate(w_sizes):
				perc_correct = [k['perc_correct'] for k in ts_perf[w_i]]
				perc_correct_sz = [k['perc_correct_sz'] for k in ts_perf[w_i]]
				plt.xlim((min(alphas) -0.01, max(alphas) +0.2  ))
				plt.ylim((min((min(perc_correct), min(perc_correct_sz))) -0.01, max((max(perc_correct), max(perc_correct_sz))) +0.01))
				plt.plot(alphas, perc_correct, 'o-r')
				plt.plot(alphas, perc_correct_sz, 's-g')
				plt.xlabel('Alpha level\n')
				plt.ylabel("Perc. Correct")
				plt.legend(labels=['Whole signal', target_str], fontsize='x-small', framealpha=0, loc='center right')
				plt.title("N=%d N_sz=%d, alphas = %s"%(ts_perf[w_i][0]['n'], ts_perf[w_i][0]['n_sz'], '  '.join([str(k) for k in alphas])))
				
				perf_filename = "%sts_%s_%d_%d_perf.png" % (save_dir, id_data, i, w_size)
				if save_dir:
					plt.savefig(perf_filename)
					plt.close()
				else:
					plt.show()
				
				# glue together the performance for this time series with the given window size and alpha levels with the rest
				# of the graphs
				all_filename = "%sts_%s_%d_%d_all.png" % (save_dir, id_data, i, w_size)
				all_filename_perf = "%sts_%s_%d_%d_all_perf.png" % (save_dir, id_data, i, w_size)
				system("convert %s %s +append %s" % (perf_filename, all_filename, all_filename_perf))

			
			# TODO put together all graphs for a given time series over all window sizes
			filenames = ' '.join(["%sts_%s_%d_%d_all_perf.png" % (save_dir, id_data, i, w_size) for w_size in w_sizes ])
			system("convert %s -append %sts_%s_%d_ALL.png" % (filenames, save_dir, id_data, i))
			


		# save histogram calculations
		results['hists'].append(hists.copy())
		hists = None
	
	
	# aggregating performance metrics
	perf = array(results['performance']['ts'])
	total =  array([ [ sum( [(lambda ts: ts['perc_correct']*ts['n'])(ts) for ts in perf[:, w_index, alpha_i]] ) 
							for alpha_i in range(len(alphas)) ] 
								for w_index in range(len(w_sizes)) ])
	n =  array([ [ sum( [ts['n'] for ts in perf[:, w_index, alpha_i]] ) 
							for alpha_i in range(len(alphas)) ] 
								for w_index in range(len(w_sizes)) ])
	print(('total', total))
	print(('n', n))
	results['performance']['total'] = total / n
	seizure =  array([ [ sum( [(lambda ts: ts['perc_correct_sz']*ts['n_sz'])(ts) for ts in perf[:, w_index, alpha_i]] ) 
							for alpha_i in range(len(alphas)) ] 
								for w_index in range(len(w_sizes)) ])
	n_sz =  array([ [ sum( [ts['n_sz'] for ts in perf[:, w_index, alpha_i]] ) 
							for alpha_i in range(len(alphas)) ] 
								for w_index in range(len(w_sizes)) ])
	print(('seizure', seizure))
	print(('n_sz', n_sz))
	results['performance']['seizure'] = seizure / n_sz
	
	
	
	"""
	# plot the performance results
	if plot:
		
		if save_dir:
			plt.savefig("%sts_%s_%d_%d_avg_std.png" % (save_dir, id_data, i, w_size))
			plt.close()
		else:
			plt.show()
	"""
	
	# TODO return the best parameters

	print(('n_bins[w_i]', n_bins[w_i])); print()
	return results




def compute_stat_perf(pvalues, target_label, grd_truth, w_size, alpha):
	"""
##
# Compute the performance for the statistical stationarity analysis
##
	"""
	# ATTENTION beware of the nan, not a number, this can contaminate the whole metric
	from numpy import average, std, sqrt, isnan, where
	from scipy.stats import chi2 
	
	from ts2cn.ts import misc

	grd_truth_w = misc.w_labels(grd_truth, w_size)

	# passed here means that we were not able to reject the null hypothesis that the windows come from the same distribution
	# it starts at the second element because the first is always -1
	passed_test = array(pvalues[1:]) > alpha
	
	# computes how much the test agrees with the labels
	n = len([p for p in pvalues if not isnan(p)]) -1
	perc_correct = sum( [# adds 1 if the two consecutives windows passed the test and the labels are equal
						# or they didn't pass but the labels changed between these windows
						1 if ( not isnan(pvalues[w_i]) ) and ( 
						passed and grd_truth_w[w_i] == grd_truth_w[w_i+1] or
						not passed and grd_truth_w[w_i] != grd_truth_w[w_i+1] )
						else 0 
						for w_i, passed in enumerate(passed_test)] ) / n
	
	# get the indices where the labels says it is a seizure
	# it starts from the second element to match with passed_test
	sz_index = where(grd_truth_w == target_label)[0][:]
	sz_index = sz_index[:-1 if sz_index[-1] == len(grd_truth_w) -1 else None] if len(sz_index) else []

	# computes how much the test agrees with the labels within the seizures
	n_sz = len(sz_index) if len(sz_index) > 0 else 1
	perc_correct_sz = sum([1 if passed_test[w_i] and grd_truth_w[w_i] == grd_truth_w[w_i+1] or
							not passed_test[w_i] and grd_truth_w[w_i] != grd_truth_w[w_i+1]
							else 0
							for w_i in sz_index]) / n_sz
	
	return {
		'perc_correct': perc_correct, 
		'perc_correct_sz': perc_correct_sz, 
		'n': n,
		'n_sz': n_sz
	}
