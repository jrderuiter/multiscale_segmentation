
import fftw3f as fftw3
import numpy as np

def fft_convolve(signal, kernel):
	s1, s2 = len(signal), len(kernel)
	size = s1 + s2 - 1
	fsize = int(2 ** np.ceil(np.log2(size)))

	paddedSignal = _zero_pad_array(signal, fsize, dtype=np.float32)
	fftSignal = fftw3.create_aligned_array(fsize//2+1, dtype=np.complex64)
	fft = fftw3.Plan(paddedSignal, fftSignal, flags=['estimate'])
	fft.execute()
	del paddedSignal

	paddedKernel = _zero_pad_array(kernel, fsize, dtype=np.float32)
	fftKernel = fftw3.create_aligned_array(fsize//2+1, dtype=np.complex64)
	fftKern = fftw3.Plan(paddedKernel, fftKernel, flags=['estimate'])
	fftKern.execute()
	del paddedKernel

	fftSignal *= (fftKernel / fsize)
	del fftKernel

	convSignal = fftw3.create_aligned_array(fsize, dtype=np.float32)
	ifft = fftw3.Plan(fftSignal, convSignal, direction=' backward', flags=['estimate'])
	ifft.execute()

	fslice = slice((s2-1)/2, (s2-1)/2+s1)
	return np.array(convSignal[fslice], dtype=np.float64)

def _zero_pad_array(arr, length, dtype=None):
	arrayLength = len(arr)
	assert arrayLength < length
 
	newArray = fftw3.create_aligned_array(length, dtype=dtype)
	newArray[:arrayLength] = arr 		# Copy over old values
	newArray[arrayLength:] = 0 			# Zero pad the rest

	return newArray