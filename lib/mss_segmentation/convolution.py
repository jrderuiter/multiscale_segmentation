'''
    Copyright (C) 2011 Julian de Ruiter
    Bioinformatics Department, Technical University of Delft, Netherlands.
	Institute for Systems Biology Seattle, Washington, USA.
	 
 	 E-Mail: j.r.deruiter@student.tudelft.nl
			 jderuite@systemsbiology.org

	Implementation of algorithm derived from MATLAB implementation
	created by T. Knijnenburg, Dutch Cancer Institute, Netherlands.
 
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.
 
    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
 
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
''' 
from scipy.ndimage import convolve1d
from scipy.signal import fftconvolve
from scipy.stats import norm

from pyfft.cl import Plan
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.elementwise import ElementwiseKernel

import sys, math, pdb

P_MIN = 1e-3
KERNEL_CUTOFF = 600

def convolve(signal, sigma, pMin=None, cutoff=None):
	if cutoff is None: cutoff = KERNEL_CUTOFF

	kernel = convolution_kernel(sigma, len(signal), pMin)	

	if len(kernel) < cutoff:
		convSignal = convolve1d(signal, kernel, mode='constant', cval=0.0)
	else:
		#convSignal = opencl_fft_1Dconvolve(signal, kernel)
		convSignal = fftconvolve(signal, kernel, mode='full')
		
		lenSignal, lenKernel = len(signal), len(kernel)
		lowerBound = (lenKernel - 1)  / 2

		convSignal = np.array(convSignal[lowerBound:lowerBound+lenSignal], copy=True, order='C')	# Copy into c-contigous array for c link compatability
	return convSignal

def opencl_fft_1Dconvolve(signal, kernel):
	ctx = cl.create_some_context(interactive=True)
	queue = cl.CommandQueue(ctx)

	s1, s2 = len(signal), len(kernel)
	size = s1 + s2 - 1
	fsize = int(2 ** np.ceil(np.log2(size)))

	paddedSig = _zero_pad_array(signal, fsize, dtype=np.complex64)
	paddedKern = _zero_pad_array(kernel, fsize, dtype=np.complex64)

	plan = Plan((fsize,1), queue=queue, dtype=np.complex64, fast_math=False)
	gpu_sig = cl_array.to_device(ctx, queue, paddedSig)
	gpu_kern = cl_array.to_device(ctx, queue, paddedKern)

	plan.execute(gpu_sig.data)
	plan.execute(gpu_kern.data)

	multKern = ElementwiseKernel(ctx, "cfloat_t *sig, cfloat_t *kern", "sig[i] = cfloat_mul(sig[i], kern[i])")
	multKern(gpu_sig, gpu_kern)

	plan.execute(gpu_sig.data, inverse=True)
	result = gpu_sig.get()

	fslice = slice((s2-1)/2, (s2-1)/2+s1)
	return np.real(result[fslice])

def _zero_pad_array(arr, length, dtype=None):
	arrayLength = len(arr)
	assert arrayLength < length
 
	newArray = np.empty(length, dtype=dtype)
	newArray[:arrayLength] = arr 		# Copy over old values
	newArray[arrayLength:] = 0 			# Zero pad the rest

	return newArray
	
def convolution_kernel(sigma, signalLength, pMin=None):
	if pMin is None: pMin = P_MIN

	normDistr = norm(loc=0, scale=sigma)
	x = round(normDistr.ppf(pMin))
	kernel = normDistr.pdf(np.arange(x, (-1 * x) + 1))

	return kernel / kernel.sum()
