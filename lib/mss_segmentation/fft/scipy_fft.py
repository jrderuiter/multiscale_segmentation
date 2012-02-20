import numpy as np

from scipy.signal import fftconvolve
from scipy.fftpack import fft, ifft

def fft_convolve(signal, kernel):
	convSignal = fftconvolve(signal, kernel, mode='full')
	s1, s2 = len(signal), len(kernel)
	fslice = slice((s2-1)/2, (s2-1)/2+s1)
	convSignal = np.array(convSignal[fslice], copy=True, order='C')	# Copy into c-contigous array
	return convSignal

def fftpack_convolve(signal, kernel):
	s1, s2 = len(signal), len(kernel)
	size = s1 + s2 - 1
	fsize = int(2 ** np.ceil(np.log2(size)))

	fftSignal = fft(signal, n=fsize)
	fftSignal *= fft(kernel, n=fsize, overwrite_x=True)
	convSignal = ifft(fftSignal, n=fsize, overwrite_x=True)

	fslice = slice((s2-1)/2, (s2-1)/2+s1)
	return np.array(convSignal[fslice], copy=True, order='C')

def _zero_pad_array(arr, length, dtype=None):
	arrayLength = len(arr)
	assert arrayLength < length
 
	newArray = np.array(length, dtype=dtype)
	newArray[:arrayLength] = arr 		# Copy over old values
	newArray[arrayLength:] = 0 			# Zero pad the rest

	return newArray