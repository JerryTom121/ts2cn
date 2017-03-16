# code to save to file and read from file the final results

# importing libraries
from json import dumps, loads

def save_to_file(path, filename, data):
	"""
--------------------------------------------
Saves the data object into the file
--------------------------------------------
path:		String.
filename:	String.
data:		Dictionary. The data to be saved. This is
			converted in json format using the dumps
			method of the json package
--------------------------------------------
Returns number of bytes written to the file 
--------------------------------------------
	"""
	file = open(path + filename, 'a')
	saved = file.write(dumps(data))
	file.close()
	return saved


def read_from_file(path, filename):
	"""
--------------------------------------------
Reads data from file
--------------------------------------------
path:		String.
filename:	String.
--------------------------------------------
Returns a dictionary of the parsed data (string) stored
into the file. This parsing is made using the loads method
from the json package
--------------------------------------------
	"""
	file = open(path + filename, 'r')
	data = loads(file.read())
	file.close()
	return data

