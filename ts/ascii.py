# code to read time series in ascii format

# importing libraries
from numpy import array



def read(file, sep='\n'):
	"""
Function to read time series values from file
--------------------------------------------
file:	Stream where from the values will be read
sep:	Separator character
--------------------------------------------
Return a list with the values of the time series read from the file
The array is a NumPy array object of int64 elements
	"""
	return array([i for i  in file.read().split(sep) if i], dtype='int64')
	
