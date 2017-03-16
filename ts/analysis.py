# code to analyse time series

# importing libraries
from numpy import array, median, mean, std


# TODO implement method that given the time serie, window size
# 	and the overlap specification, perform the full analysis 
#	This full analysis means the metric extration and aggregation
# 	Receives as parameter the path and the filename

#--------------------------------------------
# Executes a complete analysis on a given time series
#--------------------------------------------
# 
#--------------------------------------------
# Returns dictionary with the complete analysis data
# Nomenclature: the first key is the window size 
#--------------------------------------------

# TODO define the nomenclature of the return dictionary
# TODO I think this method should be implemented in the program
def complete_analysis(ts, convert_method, convert_params, indices, window_sizes):
	
	pass


# IDEA what I can do here is the analysis of redundancy, correlation between the extracted metrics
# or I create a folder to it



def window_correlations(ts, windows):
	"""
--------------------------------------------
Calculate the correlation between windows
-------------------------------------------- 
ts:			A list of numbers (the time series)
windows:	A list of window sizes. The window sizes must be ts multiples or an exception will be raised
-------------------------------------------- 
Returns a list with the mean correlation and standard deviation of the windows

Warning: This operation can be memory and CPU intensive depending on the size of the time series and the number of windows
	"""
	from scipy.stats  import pearsonr
	from numpy import array_split
	ar_ts = array(ts)
	results = []
	for w in windows:
		ar_ts = array_split(array(ts), range(w, len(ts), w))
		corr = [pearsonr(x[1], y) for x in enumerate(ar_ts) for y in ar_ts[x[0]+1:]]
		results.append((mean(corr), std(corr)))
	
	return results



def line_length(ts):
	"""
--------------------------------------------
Calculate the line lenght of a given time series
-------------------------------------------- 
ts:	A list of numbers (the time series)
-------------------------------------------- 
Return a number
	"""
	return sum([abs(ts[i+1] - ts[i]) for i in range(len(ts)-1)])



def event_w(ts, LL_med, LL_std, TF):
	"""
--------------------------------------------
Calculate the Event of a given time series (E_w)
-------------------------------------------- 
ts:		A list of numbers (the time series)
LL_med:	A number, the line length median of the time series
LL_std:	A number, the line length standard deviation of the time series
TF:		A number, the thresholding factor
-------------------------------------------- 
Return a number, 0 or 1
1 if (LL_w - LL_med) / (TF * LL_std) >= 1
0 if (LL_w - LL_med) / (TF * LL_std) < 1
"""
	return 1 if (line_length(ts) - LL_med) / (TF * LL_std) >= 1 else 0



def events(ts, w_size, LL_med, LL_std, TF):
	"""
--------------------------------------------
Calculate the Events vector of a given time series (E)
-------------------------------------------- 
ts:		A list of numbers (the time series)
w_size:	A number, the window size, in samples
LL_med:	A number, the line length median of the time series
LL_std:	A number, the line length standard deviation of the time series
TF:		A number, the thresholding factor
-------------------------------------------- 
Return a list
	"""
	return [event_w(ts[i:i+w_size], LL_med, LL_std, TF) for i in range(0, len(ts), w_size)]


def clean_events(vector, half_w_size, invert=False):
	"""
--------------------------------------------
Clean events vector (erode and dilate)
-------------------------------------------- 
vector:		A list of numbers
half_w_size:	A number, the half window size, insamples
invert:		Boolean, if True first dilate
-------------------------------------------- 
Return a list
-------------------------------------------- 
The erode block pass a sliding window through the signal, and if the 
entire window is fill with ones, then the vector element relative to 
that window is set to one, and zero otherwise.

The dilate block pass a sliding window through the signal, and if 
there is at least one non-zero element in the window, then the vector
element relativo to that window is set to one, and zero otherwise.
 
See file srep01483-s1.pdf (Supplementary information of the paper 
Automated identification of multiple seizure-related 
and interictel epileptiform event types in the EEG of mice
DOI: 10.1038 / srep01483)
--------------------------------------------
	"""
	
	size = len(vector)
	# add a half window in the beginning and in the end of vector
	# with zeros
	vector.extend([0]*(half_w_size ))
	aux_vec = [0]*(half_w_size )
	aux_vec.extend(vector)
	vector.clear()
	vector.extend(aux_vec)

	# create a vector with the size of the window at a given index
	w_sizes = []	
	w_sizes.extend( range(half_w_size +1, 2*half_w_size +1) )
	w_sizes.extend( [2*half_w_size +1]*(len(vector) - 4*half_w_size) )
	w_sizes.extend( range(2*half_w_size, half_w_size, -1) )
	
	# perform what in the article is described as erode the vector
	for k in [1,2]:
		aux_vec = [1 
			# if invert equals False, then first erode and then dilate
			if 	(k == 1 and not invert and 
				sum(vector[ i-half_w_size : i+half_w_size +1 ]) == w_sizes[i - half_w_size]) or
				(k == 2 and not invert and sum(vector[ i-half_w_size : i+half_w_size +1 ]) >= 1) or
			# if invert equal True, then first dilate and then erode
				(k == 1 and invert and sum(vector[ i-half_w_size : i+half_w_size +1 ]) >= 1) or
				(k == 2 and invert and 
				sum(vector[ i-half_w_size : i+half_w_size +1 ]) == w_sizes[i - half_w_size])
			else 0
					for i in range(half_w_size , half_w_size+size)]
		vector.clear()
		vector.extend(aux_vec)
		
	"""
	# perform what in the article is described as erode the vector
	aux_vec = [1 if sum(vector[ i-half_w_size : i+half_w_size +1 ]) == w_sizes[i - half_w_size] else 0
				for i in range(half_w_size , half_w_size+size)]
	vector.clear()
	vector.extend(aux_vec)
	
	print('ASD')
	print(vector)
	
	# perform what in the article is described as dilate the vector
	aux_vec = [1 if sum(vector[ i-half_w_size : i+half_w_size +1 ]) >= 1 else 0
				for i in range(half_w_size, half_w_size+size)]
	vector.clear()
	vector.extend(aux_vec)
	
	"""
	return vector



def dyadic_upsample(signal, even=True, pad_element=0):
	"""
--------------------------------------------
Dyadic upsampling 
-------------------------------------------- 
signal:		List, the signal to upsample
even:			Boolean, if zero-pad in the even indexes. True is the default
pad_element:	Number, the number used to pad the signal. Zero is the default
-------------------------------------------- 
Return a list
--------------------------------------------
	"""
	offset = -1 if even else 0
	[signal.insert(i, pad_element) for i in range(len(signal)+offset, offset, -1)]
	return signal



def cleanup():
	"""
--------------------------------------------
Perform the signal cleanup 
	"""
	pass


def baseline_parameters(ts, w_size, adj_tf=True):
	"""
--------------------------------------------
Calculate baseline parameters
-------------------------------------------- 
ts:				A list, the baseline time series
w_size:			A number, the window size, in samples
adj_tf:			Boolean, if the TF should be calculated
--------------------------------------------
Return a dictionary
Each one of the dictionary entries are the baseline parameters calculated
To calculate the Thresholding Factor (TF) a series of execution are performed
and its value is increased until the LL_hits is zero. 
	"""
	
	LL_med = median(ts)
	LL_std = std(ts)
	
	TF = 1
	LL_hits = 1
	if adj_tf:
		while LL_hits > 0:
			LL_hits = sum(events(ts, w_size, LL_med, LL_std, TF))
			TF += 0.5
	
	return {'LL_med': LL_med, 'LL_std': LL_std, 'TF': TF}


# TODO generalize this function, allow to pass the labels and the way to compare them
def categorize_answers(E, signal, seiz_dur, spike_amp):
	"""
--------------------------------------------
Given a event vector, categorize the seizures, spikes and abnormal
-------------------------------------------- 
E:			List. The vector of Events (a binary one)
signal:		List. The original signal used to categorize the spikes
seiz_dur:		Number. The duration, in samples, that caracterize a seizure
spike_amp:	Number. The average amplitude to caracterize a spike
-------------------------------------------- 
Return a list of strings.
Each element of the list have one of the following strings 
'sz', 'sp', 'ab', 'no' meaning seizure, spike, abnormal or 
normal respectively, refering to the category of the sample
in the respective index.
--------------------------------------------
	"""
	# at first all signal is normal
	result = ['no'] * len(E)
	
#	print('CATE')
#	print(len(E), len(result), len(signal))
#	print(E)
#	print('\n\n')

	ev_count = 0
	start = 0
	label = ''
	for i in range(len(E)):
		if E[i] == 1:
			# save the beginning of the event sequence
			start = start if ev_count != 0 else i
			ev_count += 1
#			print('E[i] == 1')
			
		elif ev_count != 0:
			
			# is this chunk a seizure ?
			if ev_count >= seiz_dur:
				label = 'sz'
			# is the average signal in this chunk greater than spike_amp ?
			elif i-start > 0 and sum(signal[start:i]) / (i-start) >= spike_amp:
				label = 'sp'
			# if it is an event but it is not a seizure or spike, it is abnormal
			elif ev_count != 0:
				label = 'ab'

			result[start:i] = [label] * ev_count
#			print(len(result), label, ev_count, start, i, len(result[start:i]))
			ev_count = 0
			
#	print('END for')
	
	# handle the degenerate case where the whole signal is an event
	# if the signal analysed is too small the 'for' above may not identify 
	# the events, so the block below is required
	
	# is this chunk a seizure ?
	if ev_count >= seiz_dur:
		label = 'sz'
	# is the average signal in this chunk greater than spike_amp ?
	elif i-start > 0 and sum(signal[start:i]) / (i-start) >= spike_amp:
		label = 'sp'
	# if it is an event but it is not a seizure or spike, it is abnormal
	elif ev_count != 0:
		label = 'ab'
	
	if ev_count:
		result[start:i] = [label] * ev_count
	
#	print(len(result), label, ev_count, start, i, len(result[start:i]))
#	print('\n')
	
	return result


def bulk_analysis(baseline, samples, w_size, wavelet, decomp_levels):
	"""
--------------------------------------------
Perform the Bulk analysis according to the paper
Automated identification of multiple seizure-related 
and interictel epileptiform event types in the EEG of mice
DOI: 10.1038 / srep01483
-------------------------------------------- 
baseline:		List, the baseline time series
samples:		List of lists, each element of the list is a 
				time series to compare with the baseline, the samples
w_size:			Number, the window size, in samples
wavelet:		String, the wavelet to be used in the 
				multiresolution decomposition. See library PyWavelets
decomp_levels:	Number, the number of levels to decompose the signal
				See library PyWavelets
--------------------------------------------
Return a list of dictionaries. Each dic entry has the LL_hits and the vector E (events)
for every time series tested (samples)
--------------------------------------------
ATTENTION 
The window size, in samples, is applied to the decomposed signal, approximation level,
so bear in mind that the desied window size as applied to the original signal should
be divided by 2**decomp_levels to the scale be preserved. E.g., if a window size of 256
samples from the original signal represents a window interval of 256 ms, then the same
window in the decomposed signal, let's say for a 3 level decomposition, should have 
256 / 2**3 = 256/8 = 32 samples
--------------------------------------------

Usage example:

import numpy as np
from datasets import bonn
from ts2cn.ts import analysis as mt
tss = bonn.load()

baseline = tss['O']
samples = tss['S']

s=samples[1:5]; levels=4; w_size=int(400/2**levels); wavelet='db4'; mt.bulk_analysis(np.hstack(baseline), s, w_size, wavelet, levels)

levels=3; w_size=int(40/2**levels); wavelet='db4'; [sum(w) for w in [[i['LL_hits'] if i['LL_hits'] == 0 else 1  for i in mt.bulk_analysis(np.hstack(baseline), s[:100], w_size, wavelet, levels)] for j in tss[0] for s in j.values()]]


#Results
#  pp=datetime.now(); levels=4; w_size=int(40/2**levels); wavelet='db4'; print(pp); [[sum(w) for w in [[i['LL_hits'] if i['LL_hits'] == 0 else 1  for i in mt.bulk_analysis(np.hstack(baseline), s[:100], w_size, wavelet, levels)] for j in tss[0] for s in j.values()]] for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp)
#[[60, 100, 46, 33, 48], [61, 100, 42, 45, 59], [45, 100, 26, 39, 53], [34, 100, 12, 34, 45], [26, 100, 10, 33, 41], [18, 100, 7, 29, 33], [11, 100, 3, 22, 31], [6, 100, 1, 26, 28], [6, 100, 1, 27, 31], [8, 100, 1, 31, 41], [6, 100, 0, 29, 35], [6, 100, 0, 33, 39], [5, 100, 0, 27, 35], [0, 100, 0, 6, 21], [0, 99, 0, 2, 15]]
#0:19:06.277569

# pp=datetime.now(); levels=4; w_size=int(40/2**levels); wavelet='db4'; all_samples=np.vstack([s for j in tss[0] for s in j.values()]); print(pp); oo = lambda l, s : [sum(l[i : i+s]) for i in range(0, len(l), s)]; [oo([i['LL_hits'] if i['LL_hits'] == 0 else 1 for i in mt.bulk_analysis(np.hstack(baseline), all_samples, w_size, wavelet, levels)], int(len(all_samples)/5)) for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp);
#2015-07-23 16:54:14.897655
#[[60, 100, 46, 33, 48], [61, 100, 42, 45, 59], [45, 100, 26, 39, 53], [34, 100, 12, 34, 45], [26, 100, 10, 33, 41], [18, 100, 7, 29, 33], [11, 100, 3, 22, 31], [6, 100, 1, 26, 28], [6, 100, 1, 27, 31], [8, 100, 1, 31, 41], [6, 100, 0, 29, 35], [6, 100, 0, 33, 39], [5, 100, 0, 27, 35], [0, 100, 0, 6, 21], [0, 99, 0, 2, 15]]
#0:03:35.931310



# pp=datetime.now(); levels=3; w_size=int(40/2**levels); wavelet='db4'; print(pp); [[sum(w) for w in [[i['LL_hits'] if i['LL_hits'] == 0 else 1  for i in mt.bulk_analysis(np.hstack(baseline), s[:100], w_size, wavelet, levels)] for j in tss[0] for s in j.values()]] for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp)
#[[1, 80, 0, 0, 3], [1, 79, 0, 0, 3], [0, 80, 0, 0, 2], [3, 82, 0, 0, 2], [1, 78, 0, 0, 2], [1, 76, 0, 0, 2], [2, 77, 0, 0, 2], [0, 74, 0, 0, 2], [0, 79, 0, 0, 2], [0, 80, 0, 0, 2], [2, 85, 0, 0, 3], [3, 84, 0, 0, 3], [1, 80, 0, 0, 2], [2, 85, 0, 0, 2], [1, 86, 0, 0, 3]]
#1:17:53.930201

# pp=datetime.now(); levels=3; w_size=int(40/2**levels); wavelet='db4'; all_samples=np.vstack([s for j in tss[0] for s in j.values()]); print(pp); oo = lambda l, s : [sum(l[i : i+s]) for i in range(0, len(l), s)]; [oo([i['LL_hits'] if i['LL_hits'] == 0 else 1 for i in mt.bulk_analysis(np.hstack(baseline), all_samples, w_size, wavelet, levels)], int(len(all_samples)/5)) for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp);
#2015-07-23 16:57:50.850142
#[[1, 80, 0, 0, 3], [1, 79, 0, 0, 3], [0, 80, 0, 0, 2], [3, 82, 0, 0, 2], [1, 78, 0, 0, 2], [1, 76, 0, 0, 2], [2, 77, 0, 0, 2], [0, 74, 0, 0, 2], [0, 79, 0, 0, 2], [0, 80, 0, 0, 2], [2, 85, 0, 0, 3], [3, 84, 0, 0, 3], [1, 80, 0, 0, 2], [2, 85, 0, 0, 2], [1, 86, 0, 0, 3]]
#0:14:44.256077



# pp=datetime.now(); levels=2; w_size=int(40/2**levels); wavelet='db4'; print(pp); [[sum(w) for w in [[i['LL_hits'] if i['LL_hits'] == 0 else 1  for i in mt.bulk_analysis(np.hstack(baseline), s[:100], w_size, wavelet, levels)] for j in tss[0] for s in j.values()]] for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp)
#2015-07-22 15:52:17.661901
[[0, 78, 0, 0, 2], [2, 75, 0, 0, 2], [1, 71, 0, 0, 0], [0, 72, 0, 0, 0], [1, 72, 0, 0, 1], [1, 73, 0, 0, 0], [0, 71, 0, 0, 0], [0, 72, 0, 0, 0], [1, 71, 0, 0, 0], [0, 71, 0, 0, 0], [0, 70, 0, 0, 0], [3, 74, 0, 0, 0], [1, 75, 0, 0, 0], [1, 73, 0, 0, 0], [2, 74, 0, 0, 0]]
#2:12:02.357042

# pp=datetime.now(); levels=2; w_size=int(40/2**levels); wavelet='db4'; all_samples=np.vstack([s for j in tss[0] for s in j.values()]); print(pp); oo = lambda l, s : [sum(l[i : i+s]) for i in range(0, len(l), s)]; [oo([i['LL_hits'] if i['LL_hits'] == 0 else 1 for i in mt.bulk_analysis(np.hstack(baseline), all_samples, w_size, wavelet, levels)], int(len(all_samples)/5)) for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp);
#2015-07-23 17:12:35.156734
#[[0, 78, 0, 0, 2], [2, 75, 0, 0, 2], [1, 71, 0, 0, 0], [0, 72, 0, 0, 0], [1, 72, 0, 0, 1], [1, 73, 0, 0, 0], [0, 71, 0, 0, 0], [0, 72, 0, 0, 0], [1, 71, 0, 0, 0], [0, 71, 0, 0, 0], [0, 70, 0, 0, 0], [3, 74, 0, 0, 0], [1, 75, 0, 0, 0], [1, 73, 0, 0, 0], [2, 74, 0, 0, 0]]
#0:27:24.785785
	"""
	from ts2cn.pre_processing import multiresolution as multires
	
	# decompose the signal, and only return the approximation level
	baseline_decomp = multires.decompose(baseline, decomp_levels, wavelet, mode='sp1')[0]
	samples_decomp = [multires.decompose(i, decomp_levels, wavelet, mode='sp1')[0] for i in samples]
	
	# extract the baseline parameters
	params = baseline_parameters(ts=baseline_decomp, w_size=w_size)
	
	# test every sample, events vector for every time series
	E = [events(ts, w_size, params['LL_med'], params['LL_std'], params['TF'])
			for ts in samples_decomp]
	
	# return the LL_hits and the events vector for every time series
	return [{
		'LL_hits': sum(i),
		'E': i
	} for i in E]


# TODO convert all lists to np.array

def automated_categorization(baseline, samples, labels_vector, w_size, wavelet, decomp_levels, 
		seiz_dur, spike_amp, clean_event_window=2, TF_range=[i/2 for i in range(1, 20)], 
		plot_graphs='save', plot_labels='', signal_ylim=None, directory='/tmp/', silent=True, 
		labels_str=['no', 'sz', 'sp', 'ab'], labels_titles=['Normal', 'Seizure', 'Spikes', 'Abnormal'], 
		labels_pos=-600,
		pred_style=['red', 'green', 'blue', 'magenta'], grd_style=['red', 'green', 'blue', 'magenta'],
		aggregate_labels=[['no'], ['sz', 'sp', 'ab']], aggregate_titles=['Normal', 'Events'], 
		identification=''):
	"""
--------------------------------------------
Perform the automated categorization of EEG events according to the paper
Automated identification of multiple seizure-related 
and interictel epileptiform event types in the EEG of mice
DOI: 10.1038 / srep01483
-------------------------------------------- 
baseline:			List, the baseline time series
samples:			List of lists, each element of the list is a 
					time series to compare with the baseline, the samples
labels_vector:		List of lists, the label for each window of each sample (Time Series)
					For each sample with N elements, it should have int(N/w_size) values
w_size:				Number, the window size, in samples
wavelet:			String, the wavelet to be used in the 
					multiresolution decomposition. See library PyWavelets
decomp_levels:		Number, the number of levels to decompose the signal
					See library PyWaveletsa
seiz_dur:     		Number. The duration, in number of windows (w_size), 
					that caracterize a seizure
spike_amp:    		Number. The average amplitude to caracterize a spike
clean_event_window:	Number, the half window size, in samples, for eroding 
					and dilating the Events vector. See function clean_event_window.
TF_range:			List, the values used as TF (Thresholding Factor)
plot_graphs:		String. 'save' to save to file the plots generated, 'plot' to display the generated graphs
					None to not display or save graphs
plot_labels:		String or List. A identification string for the plots
signal_ylim:		Tuple or None. Set the y limit for signal plot
directory:			String. The directory where to save the generated plots if requested
silent:				Boolean. Used to debug, default True
labels_str:			List of string. The string representing each class (label) present in the data
labels_titles:		
label_pos:			Number. The x position of the text "Predicted" and "Ground Truth" of the signal plot
pred_style:			
grd_style:			
aggregate_labels:	
aggregate_titles:	
identification:		String. Just a string to identify the execution
-------------------------------------------- 
 
Return a list of dictionaries. Each dic entry has the LL_hits and the vector E (events)
for every time series tested (samples) so bear in mind that the desied window size 
as applied to the original signal should
be divided by 2**decomp_levels to the scale be preserved. E.g., if a window size of 256
samples from the original signal represents a window interval of 256 ms, then the same
window in the decomposed signal, let's say for a 3 level decomposition, should have 
256 / 2**3 = 256/8 = 32 samples
--------------------------------------------
ATTENTION
Due to effects of the wavelet decomposition and the window size not necessarily being a multiple
of the signal's length the final vector has a size different from the original signal.
To address this effect the final vector is truncated to the size of the labels_vector.
--------------------------------------------

Usage example:

import numpy as np
import imp
from matplotlib import pyplot as plt
from ts2cn.ts import analysis as mt
imp.reload(plt)
from datetime import datetime
from datasets import bonn
from ts2cn.ts import analysis as mt
tss = bonn.load()

baseline = tss['O']
samples = tss['S']

pp=datetime.now(); s=samples[1:5]; levels=4; w_size=int(400/2**levels); wavelet='db4'; labels = [['sz'] * len(s[0]/w_size)] * len(s); seiz_dur=1000; spike_amp=1000; clean_event_window=8; TFs=[1,5,10]; print(pp); mt.automated_categorization(np.hstack(baseline), s, labels, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window, TFs); print(datetime.now() - pp);




# articles parameters
seiz_dur = int(5*173.61)					# 5 seconds
spike_amp = 250 							# 250 uV
clean_event_window = int(0.040 * 173.61)	# 40 miliseconds
TFs = [i/2 for i in range(1, 21)]


all_samples=np.vstack([s for j in tss[0] for s in j.values()]); 
labels1 = [['no'] * 100 , ['sz'] * 100, ['no'] * 100, ['no'] * 100, ['sp'] * 100]
labels2 = [['no'] * 100 , ['sz'] * 100, ['no'] * 100, ['sp'] * 100, ['sp'] * 100]
labels3 = [['no'] * 100 , ['sz'] * 100, ['no'] * 100, ['ab'] * 100, ['sp'] * 100]
labels4 = [['no'] * 100 , ['sz'] * 100, ['no'] * 100, ['no'] * 100, ['no'] * 100]

pp=datetime.now(); levels=3; w_size=int(40/2**levels); wavelet='db4'; print(pp); [[i for i in mt.automated_categorization(np.hstack(baseline), all_samples, labels1, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window, TFs)], int(len(all_samples)/5) for w_size in [10, 20, 30, 40, 50, 60, 80, 120, 160, 200, 240, 280, 320, 360, 400]]; print(datetime.now() - pp);



# teste de seizure
data = tss['S'][:2]
data = tss['F'][:2]

pp=datetime.now(); levels=3; w_size=int(40/2**levels); labels = [['sz']*int(len(data[0])/w_size)] * len(data); wavelet='db4'; print(pp); oo = mt.automated_categorization(np.hstack(baseline), data, labels, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window, TFs[:3], silent=True)


data=tss['S'][:8]; FS=173.61; imp.reload(mt); pp=datetime.now(); levels=3; w_size=int(0.240*FS/2**levels); seiz_dur=int(5*173.61/w_size); labels = [['sz']*int(len(data[0])/w_size)] * len(data); wavelet='db4'; print(pp); oo = mt.automated_categorization(np.hstack(baseline), data, labels, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window=1, TF_range=TFs[:20], silent=True)



data=list(np.vstack((tss['O'][:100], tss['S'][:100]))); FS=173.61; imp.reload(mt); pp=datetime.now(); levels=3; w_size=int(0.240*FS/2**levels); seiz_dur=int(5*173.61/w_size); labels = list(np.vstack(([['no']*int(len(data[0])/w_size)] * int(len(data)/2), [['sz']*int(len(data[0])/w_size)] * int(len(data)/2)))); wavelet='db4'; print(pp); oo = mt.automated_categorization(np.hstack(baseline), data, labels, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window=1, TF_range=TFs[:20], silent=True)



data=list(np.vstack((tss['O'][:100], tss['S'][:100]))); FS=173.61; imp.reload(mt); pp=datetime.now(); levels=3; w_size=int(0.240*FS/2**levels); seiz_dur=int(5*173.61/w_size); labels = list(np.vstack(([['no']*int(len(data[0])/w_size)] * int(len(data)/2), [['sz']*int(len(data[0])/w_size)] * int(len(data)/2)))); wavelet='db4'; print(pp); oo = mt.automated_categorization(np.hstack(baseline), data, labels, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window=1, TF_range=TFs[:20], silent=True, plot_graphs='save', directory='/tmp/automated/', identification='Teste de execucao')
	"""
	
	from ts2cn.pre_processing import multiresolution as multires
	from ts2cn.plot import graphics
	from sklearn.metrics import accuracy_score, f1_score, precision_score,\
								recall_score, roc_curve, roc_auc_score
	from numpy import hstack, ndarray, array, argmax
	from datetime import datetime
	import imp, os

	imp.reload(multires)
	
	# Prevent the print function of printing
	global print
	
	p_function = print
	def silent_print(*args, **kwargs):
		pass
	if silent:
		print = silent_print

	print('Começando ' + str(datetime.now()))
	
	# Decompose the signal
	baseline_decomp = multires.decompose(baseline, levels=decomp_levels, wavelet=wavelet, mode='sp1')[0]
	samples_decomp = [multires.decompose(array(i), levels=decomp_levels, wavelet=wavelet, mode='sp1')[0] for i in samples]

	print('Decomposição concluída ' + str(datetime.now()))
	
	# Extract baseline parameters
	params = baseline_parameters(ts=baseline_decomp, w_size=w_size, adj_tf=False)

	print(params)
	print('Parametros de baseline encontrados ' + str(datetime.now()))
	
	# Iterate over TF_range (Thresholding Factor values)
	# Apply the clean_events function, to erode and dilate
	
	# For every TF in TF_range create a events vector (E), so the first dimension of the TF_E
	# vector is related to the TF value used to calculate the events vector (E)
	TF_E = [[clean_events( events(ts, w_size, params['LL_med'], params['LL_std'], TF),
				clean_event_window) for ts in samples_decomp] for TF in TF_range]
	
#	print('Sinal sem limpar '+ str(datetime.now()))
#	print([[events(ts, w_size, params['LL_med'], params['LL_std'], TF)
#				for ts in samples_decomp] for TF in TF_range] )

#	print(TF_E)
#	print('Calculando e limpando eventos ' + str(datetime.now()))
	
	# Upsample the signal, the number of decomposition levels padding with zero
	# and in the even indexes
	for i in range(len(TF_E)):
		for j in range(len(TF_E[i])):
			
			for k in range(decomp_levels):
				TF_E[i][j] = dyadic_upsample(TF_E[i][j], even=True, pad_element=0)
		
	
#	print(TF_E)
#	print('Upsampling complete ' + str(datetime.now()))
	
	# dilate and then erode after the upsampling as described in the paper, using a 
	# window (2**decomp_levels-1) the size of the previous one because of the upsampling
	# and the "-1" because it is a half-window
	TF_E = [[clean_events(E, int(clean_event_window*2**(decomp_levels-1)), invert=True) for E in TF] for TF in TF_E]

	print(TF_E)
	print('\n\nTAMANHOS')
	print([[len(E) for E in TF] for TF in TF_E])
	
	# Create a answers vector, has the same dimension as the vector TF_E
	# Also truncate the vector size to be the same as labels_vector
	ans = [[categorize_answers(TF[ev], samples[ev], seiz_dur, spike_amp)[:len(labels_vector[0])] 
			for ev in range(len(TF))] for TF in TF_E]
	
	print([[(ev, len(TF[ev]), len(samples[ev])) 
			for ev in range(len(TF))] for TF in TF_E])
	print('\n\nTAMANHOS')
	print([[len(E) for E in TF] for TF in ans])
	print('Categorizando respostas ' + str(datetime.now()))

	print(len(labels_vector))
	print(len(labels_vector[0]))
	
	
	if silent:
		print = p_function

	#return ans
	
	# compute the performance metrics, accuracy, recall, precision and F-Score
	# TODO check if it is better to analyse all together like a big vector concatenating all
	# samples in a single vector or the way it is, analysing per sample
	# If consider as a big vector, do this transformation when creating the vector ans above
	# No it is better to return the vector ans as it is and a vector with the metrics per TF
	# i.e., for each TF executed merge all windows and compute the measures below as a whole 
	# for this given TF
	ground_truth = hstack(labels_vector)
	TFs = array([hstack(TF) for TF in ans])
	
	# DEBUG REMOVE LATER
	print('BEFORE ERROR')
	print((len(ground_truth), ground_truth))
	print((len(TFs), len(TFs[0]), TFs))
	print((len(labels_str), labels_str))
	
	results = [{
		'f1': f1_score(ground_truth, TF, labels=labels_str, average=None),
		'accuracy': accuracy_score(ground_truth, TF),
		'recall': recall_score(ground_truth, TF, labels=labels_str, average=None),
		'precision': precision_score(ground_truth, TF, labels=labels_str, average=None)
	} for TF in TFs]
	
	
	# aggregate labels and compute the metrics in this case
	aggregate_TFs = TFs.copy()
	aggregate_grd = ground_truth.copy()
	# replace the orignal labels by its position in the array aggregate_labels
	for new_label in range(len(aggregate_labels)):
		for old_label in aggregate_labels[new_label]:
			
			for TF in range(len(aggregate_TFs)):
				aggregate_TFs[TF][ aggregate_TFs[TF] == old_label ] = new_label
			aggregate_grd[ aggregate_grd == old_label ] = new_label
	
	aggregate_labels_str = [str(i) for i in range(len(aggregate_labels))]
	aggregate_results = [{
		'f1': f1_score(aggregate_grd, TF, labels=aggregate_labels_str, average=None),
		'accuracy': accuracy_score(aggregate_grd, TF),
		'recall': recall_score(aggregate_grd, TF, labels=aggregate_labels_str, average=None),
		'precision': precision_score(aggregate_grd, TF, labels=aggregate_labels_str, average=None)
	} for TF in aggregate_TFs]
	
	print('Computation complete')


	################################################
	################################################
	# plot graphs
	################################################
	################################################
	if plot_graphs:
		imp.reload(graphics)
		from matplotlib import pyplot as plt
		
		execution_date = datetime.now()
		
		# saving execution context 
		if not os.path.exists(directory):
		    os.makedirs(directory)

		file = open(directory + 'execution.txt', 'a')
		string = "date: %s\n" % execution_date
		string += 'identification = ' + identification
		string += "\n(baseline[:20], samples[:20], labels_vector[:20], w_size, wavelet, decomp_levels, \
			seiz_dur, spike_amp, clean_event_window, TF_range, \
			plot_graphs, directory, silent, labels_str)\n"
		context = (str(baseline)[:300], str(samples)[:300], str(labels_vector)[:300], w_size, wavelet, decomp_levels, 
			seiz_dur, spike_amp, clean_event_window, TF_range, 
			plot_graphs, directory, silent, labels_str)
		string += '\n'.join([str(i) for i in context])
		string += "\n#################\n"
		file.write(string)
		file.close()
		
		
		print('Generating plots - %s' % execution_date)
		
		
		
		################################################
		# plot all scores for every threshold factor
		################################################
		scores_dir='scores/'; scores_figsize=(8, 6);
		if not os.path.exists(directory + scores_dir):
		    os.makedirs(directory + scores_dir)
		scores_labels = ['Accuracy', 'Precision', 'Recall', 'F1']
		# TODO allow customization of colors
		scores_style = [['r', 'o', '-'], ['b', '*', '-'], ['g', 'v', '-'], ['y', '^', '-']]
		
		for index in range(len(labels_titles)):
			# get a array with the score for the respective metric (accuracy, precision, recall, or f1)
			# this is needed because the accuracy is a scalar while the other metrics are arrays
			scores = [tuple(j[index] if type(j) == type(ndarray([1])) else j 
						for j in i.values()) for i in results]; 
			
			# plot the scores
			graphics.plot_scores(scores, x=TF_range, labels=scores_labels, x_label='Threshold Values', 
				y_label='Score', style=scores_style, leg_loc='best', filename=labels_titles[index]+'.png', 
				directory=directory + scores_dir, figsize=scores_figsize, display=plot_graphs);

		for index in range(len(aggregate_labels)):
			# aggregate scores 
			scores_aggregate = [tuple(j[index] if type(j) == type(ndarray([1])) else j 
						for j in i.values()) for i in aggregate_results]; 
			
			# plot aggregate scores
			graphics.plot_scores(scores_aggregate, x=TF_range, labels=scores_labels, x_label='Threshold Values', 
				y_label='Score', style=scores_style, leg_loc='best', filename='aggregate_' + aggregate_titles[index]+'.png', 
				directory=directory + scores_dir, figsize=scores_figsize, display=plot_graphs);
		
		# save scores to file
		file = open(directory + scores_dir + 'scores.txt', 'a')
		string = "date: %s\n" % execution_date
		string += "TF: %s\n" % TF_range
		string += "labels: %s\n" % labels_str
		string += "scores: \n%s\n\n" % '\n'.join([str(i) for i in results])
		string += "aggregate_labels: %s\n" % aggregate_labels_str
		string += "aggregate_scores: \n%s\n" % '\n'.join([str(i) for i in scores_aggregate])
		string += "#################\n"
		file.write(string)
		file.close()
		
		################################################
		# plot confusion matrix, also save to file
		################################################
		from sklearn.metrics import confusion_matrix
		
		conf_mat = [confusion_matrix(ground_truth, TF, labels_str) for TF in TFs]
		conf_mat_aggregate = [confusion_matrix(aggregate_grd, TF, aggregate_labels_str) for TF in aggregate_TFs]
		
		conf_mat_dir = 'conf_mat/'
		if not os.path.exists(directory + conf_mat_dir):
		    os.makedirs(directory + conf_mat_dir)
#		pred_style = ['red', 'green', 'blue', 'magenta']
#		grd_style = ['red', 'green', 'blue', 'magenta']
		grd_text = 'Ground Truth'
		pred_text = 'Predicted'
		conf_mat_fig_title = 'TF = %.1f' # the title of the figure
		conf_mat_matrix_title = 'Confusion Matrix'
		conf_mat_plot_title = '' # the title of the plot, the one at the left
		conf_mat_ylabel = 'Percentage of Ground Truth'
		conf_mat_figsize = (8, 6)
		for i in range(len(conf_mat)):
			graphics.plot_conf_mat(conf_mat[i], labels=labels_str, pred_style=pred_style, grd_style=grd_style, 
					grd_text=grd_text, pred_text=pred_text, bar_width=0.4, bar_offset=0.8, title=conf_mat_fig_title % TF_range[i], 
					display=plot_graphs, filename="conf_mat_TF_%s.png" % TF_range[i], directory=directory + conf_mat_dir, 
					figsize=conf_mat_figsize, wspace=0.2, hspace=0.2, top=0.95, bottom=0.2, right=0.95, left=0.1, 
					mat_title=conf_mat_matrix_title, plot_title=conf_mat_plot_title, ylabel_plot=conf_mat_ylabel)
			
			graphics.plot_conf_mat(conf_mat_aggregate[i], labels=['no', 'ev'], pred_style=pred_style[:len(aggregate_labels_str)], 
					grd_style=grd_style[:len(aggregate_labels_str)], grd_text=grd_text, pred_text=pred_text, bar_width=0.2, 
					bar_offset=0.5, title=conf_mat_fig_title % TF_range[i], display=plot_graphs, filename="aggregate_conf_mat_TF_%s.png" % TF_range[i], 
					directory=directory + conf_mat_dir, figsize=conf_mat_figsize, wspace=0.2, hspace=0.2, top=0.95, bottom=0.2, 
					right=0.95, left=0.1, mat_title=conf_mat_matrix_title, plot_title=conf_mat_plot_title, ylabel_plot=conf_mat_ylabel)
			
		# save confusion matrix to file
		file = open(directory + conf_mat_dir + 'conf_mat.txt', 'a')
		string = "date: %s\n" % execution_date
		string += "TF: %s\n" % TF_range
		string += "conf_mat: \n%s\n" % '\n'.join([str({i: conf_mat[i]}) for i in range(len(conf_mat))])
		string += "aggregate_conf_mat: \n%s\n" % '\n'.join([str({i: conf_mat_aggregate[i]}) for i in range(len(conf_mat_aggregate))])
		string += "#################\n"
		file.write(string)
		file.close()
			
		################################################
		# print to file and screen statistics, odds ratio and Chi^2
		################################################
		# ATTENTION only calculate this statistics for the aggregate values (Normal and events)
		from scipy.stats import chi2_contingency
		from scipy.stats import fisher_exact
		# computing Chi^2
		chi2 = []
		odds = []
		# ATTENTION there some cases where only one label is classified, then the chi2_contingency got an exception
		# The same happens when generating the inviduals ROC curves, so the index there is based on the chi2 vector
		for mat in conf_mat_aggregate:
			try:
				chi2.append(chi2_contingency(mat)[:3])
				odds.append(fisher_exact(mat))
			except ValueError as e:
				print(e)

		 
		stats_dir = 'statistics/'
		if not os.path.exists(directory + stats_dir):
		            os.makedirs(directory + stats_dir)
		# save statistics to file
		file = open(directory + stats_dir + 'statistics.txt', 'a')
		string = "date: %s\n" % execution_date
		string += "TF: %s\n" % TF_range
		string += "Chi^2 (X^2, p-value, degree of freedom) \n%s\n" % '\n'.join([str({i: chi2[i]}) for i in range(len(chi2))])
		string += "odds_ratio: \n%s\n" % '\n'.join([str({i: odds[i]}) for i in range(len(odds))])
		string += "#################\n"
		file.write(string)
		file.close()
		
		print('Statistics')
		print(string)
		
		################################################
		# plot ROC curves, all of them, aggregate and from classifier
		################################################
		roc_dir = 'roc/'
		curves_dir = 'curves/'
		if not os.path.exists(directory + roc_dir):
			os.makedirs(directory + roc_dir)
		if not os.path.exists(directory + roc_dir + curves_dir):
			os.makedirs(directory + roc_dir + curves_dir)
		roc_figsize = (8, 6)
		# aggregate with complement
		graphics.plot_roc_curve(predicted=ans, ground_truth=[labels_vector] * len(ans), labels=labels_str, aggregate=True, 
					compl={'no':True, 'sz':False, 'sp':False, 'ab':False}, titles=['All Events'] + labels_titles[1:], 
					figsize=roc_figsize, hspace=0.4, wspace=0.3, leg_fontsize='small', suptitle="Roc Curves", 
					display=plot_graphs, filename="aggregate_roc.png", directory=directory + roc_dir)

		# plot for each TF, with complement
		[graphics.plot_roc_curve(predicted=ans[i], ground_truth=labels_vector, labels=labels_str, aggregate=False, 
					compl={'no':True, 'sz':False, 'sp':False, 'ab':False}, titles=['All Events'] + labels_titles[1:], 
					figsize=roc_figsize, hspace=0.4, wspace=0.3, leg_fontsize='small', 
					suptitle="Roc Curve TF = %.1f" % TF_range[i], display=plot_graphs, filename="roc_TF_%s.png" % TF_range[i], 
					directory=directory + roc_dir + curves_dir) 
					for i in range(len(chi2))]
		
		
		
		################################################
		# plot signal with predicted and ground truth values
		################################################
		
		signal_dir = 'signal/'
		if not os.path.exists(directory + signal_dir):
			os.makedirs(directory + signal_dir)
		signal_figsize = (15, 5)
		# TODO customize the colors
		signal_labels_style = {'no':'darkblue', 'ab':'yellow', 'sz':'orange', 'sp':'maroon'} #grd_style #['red', 'blue', 'green', 'magenta']
		signal_labels_dic = {'ab': 'Abnormal', 'no': 'Normal', 'sp':'Spike', 'sz':'Seizure'}
		sig_lbls_stl = [signal_labels_style[i] for i in signal_labels_dic.keys()]
		signal_title = "Second by Second Analysis - TF = %.1f"
		
		if identification:
			signal_title = identification +' - '+ signal_title
		
		min_bar_height = 2
		
		#ymax = max([max(sample) for sample in samples])
		#ymin = min([min(sample) for sample in samples])
		
		# only plot for five executions, two above and two below the best chi squared measure, centered 
		# in the one with the best chi^2 measure (highest)
		# ATTENTION, this is just a heuristic
		best_chi2 = 0
		if chi2:
			best_chi2 = argmax([i[0] for i in chi2])
		
		imp.reload(plt)
		imp.reload(graphics)
		
		# The range goes from the max(best_chi2 -2) to min(best_chi2 +3) because
		# the first can be negative and the later can be greater than len(TF_range)
		for TF in range(max(best_chi2 -2, 0), min(best_chi2 + 3, len(TF_range))):
			if not os.path.exists(directory + signal_dir + 'TF_' + str(TF_range[TF])):
				os.makedirs(directory + signal_dir + 'TF_' + str(TF_range[TF]))
			for i in range(len(samples)):
				# sets the bar height to 5% of the range of the signal
				bar_height = max(int( ((max(samples[i]) - min(samples[i]))/100)*5 ), min_bar_height)
				
				signal_label = ' Signal'
				plot_label_str = ''
				if type(plot_labels) == type(''):
					plot_label_str = plot_labels
				else:
					plot_label_str = plot_labels[i]
				signal_label = plot_label_str + signal_label
				
				graphics.plot_sig_pred_grd(samples[i], 'Samples', signal_label, 'c-', labels_vector[i], 'Ground Truth',  ans[TF][i], 'Predicted', 
					signal_labels_dic, sig_lbls_stl, labels_pos=labels_pos, window=w_size, bar_height=bar_height, v_offset=0, ylim=signal_ylim, 
					title=signal_title % TF_range[TF], display=plot_graphs, directory=directory + signal_dir + 'TF_' + str(TF_range[TF]) +'/', 
					filename='sample_%s_%s.png' % (i, plot_label_str), figsize=signal_figsize)
			
	
	return {'metrics':results, 'metrics_aggregate':aggregate_results, 'events': ans}


