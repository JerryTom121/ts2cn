# Perform a multiresolution decomposition of the signal (time series)
# Based on the library PyWavelets http://pywavelets.readthedocs.org/en/latest/index.html

import numpy as np
from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt
import pywt

# TODO allow to save to file
# TODO move this function to the filter file
def print_freq_resp(ftr):
	"""
----------------------------------------------------------------
Plot the frequency response of a given filter
---------------------------------------------------------------- 
ftr:	The filter object, same returned from butter function
----------------------------------------------------------------
	"""
	w, h = freqz(ftr)
	fig = plt.figure()
	plt.title('Digital filter frequency response')
	ax1 = fig.add_subplot(111)
	plt.plot(w, 20 * np.log10(abs(h)), 'b')
	plt.ylabel('Amplitude [dB]', color='b')
	plt.xlabel('Frequency [rad/sample]')
	ax2 = ax1.twinx()
	angles = np.unwrap(np.angle(h))
	plt.plot(w, angles, 'g')
	plt.ylabel('Angle (radians)', color='g')
	plt.grid()
	plt.axis('tight')
	plt.show()	



def decompose(signal, levels, wavelet, mode='sp1'):
	"""
----------------------------------------------------------------
Perform a multiresolution decomposition
---------------------------------------------------------------- 
signal:	The signal to decompose
levels:	The number of levels to decompose
wavelet:	The wavelet name. See pyWavelets documentation
			http://pywavelets.readthedocs.org/en/latest/ref/wavelets.html#wavelet-object
mode:		the mode of the convolution, see pyWavelets documentation
			http://pywavelets.readthedocs.org/en/latest/ref/signal-extension-modes.html#ref-modes
			options:
 			zpd - zero-padding - signal is extended by adding zero samples
 			cpd - constant-padding - border values are replicated
 			sym - symmetric-padding - signal is extended by mirroring samples
 			ppd - periodic-padding - signal is treated as a periodic one
 			sp1 - smooth-padding - signal is extended according to the first derivatives calculated on the edges (straight line)
---------------------------------------------------------------- 
Return a list of lists, where each element of the first list is the decomposed coefficients.
The order of the coefficients is backwards, the last approximation level first, then the
last detail level and so on, the last element of the list is the coefficients for the first detail.
	"""
	sig = np.array(signal)
	
	return pywt.wavedec(sig, wavelet, mode=mode, level=levels)


def reconstruct(decomp, wavelet, mode='sp1'):
	"""
----------------------------------------------------------------
Reconstruct the signal from the decomposed components
---------------------------------------------------------------- 
decomp:	A list with the decomposed components. See decompose function
wavelet:	The wavelet name. See pyWavelets documentation
			http://pywavelets.readthedocs.org/en/latest/ref/wavelets.html#wavelet-object
mode:		the mode of the convolution, see pyWavelets documentation
			http://pywavelets.readthedocs.org/en/latest/ref/signal-extension-modes.html#ref-modes
			options:
 			zpd - zero-padding - signal is extended by adding zero samples
 			cpd - constant-padding - border values are replicated
 			sym - symmetric-padding - signal is extended by mirroring samples
 			ppd - periodic-padding - signal is treated as a periodic one
 			sp1 - smooth-padding - signal is extended according to the first derivatives calculated on the edges (straight line)
---------------------------------------------------------------- 
Return a list with the reconstructed signal
	"""
	return np.array(pywt.waverec(decomp, wavelet, mode))


# Usage example
# imp.reload(multiresolution); wavelet = 'db4'; levels = 4; y4 = multiresolution.decompose(data, levels, wavelet); multiresolution.print_decomposition(data, y4, fs); data4 = multiresolution.reconstruct(y4, wavelet); plt.plot(data, label='Orig'); plt.plot(data4, label='Recon'); plt.legend(); plt.show(); print('Erro abs de reconstrucao ' + str(sum([abs(data[i] - data4[i]) for i in range(len(data))])));



def print_decomposition(orig_sig, decom_sig, fs):
	# plot the original signal
	levels = len(decom_sig)
	plt.subplot(levels+1, 1, 1)
	plt.plot(orig_sig, 'b')
	plt.title("Original Signal")
	plt.grid()
	
	# plotting the detail levels
	for l in range(levels-1, -1, -1):
		plt.subplot(levels+1, 1, levels-1 - l+2)
		plt.plot(decom_sig[l], 'b')
		plt.title("Detail level D"+ str(levels - l))
		plt.grid()
		plt.legend()
		
	# plotting the last approximation level
	plt.subplot(levels+1, 1, levels+1)
	plt.plot(decom_sig[0], 'b')
	plt.title("Approx. Level A"+ str(levels-1))
	plt.grid()
	plt.legend()

	# space between subplots
	plt.subplots_adjust(hspace=1)
	
	plt.show()

