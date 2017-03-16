# Miscelanea code to time series analysis

# importing libraries
from numpy import array





def ts_windowed_indices(ts_size, window_size=100, overlap=0):
	"""
##
# Generate the index used in the convertion of
# time series through the windowing method
##
# ts_size:		The length of the time series
# window_size:	The window size
# overlap:		The number of observations the current window overlap
#				with the previous window
##
# return a list of index where the sub-series start 
##
	"""
	last_index = ts_size - window_size 
	return [window_size*i - i*overlap for i in range(int(last_index/(window_size-overlap))+1)]


def w_labels(labels, w_size):
	"""
##
# Compute the labels per window
##
# labels:	List. The list of labels, it should have one label instance per item on the time series
# w_size:	Number. The window size
##
# Return an array with the most frequent, mode, for each window
##
	"""
	return array([ max( set(labels[ (k*w_size):((k+1)*w_size) ]), key=list(labels[ (k*w_size):((k+1)*w_size) ]).count ) 
					for k in range(int(len(labels)/w_size)) ])
