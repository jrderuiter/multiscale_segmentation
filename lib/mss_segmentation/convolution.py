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
from numpy import arange, array
from scipy.ndimage import convolve1d
from scipy.stats import norm

from fft.fftw import fft_convolve as fftw_convolve
from fft.scipy_fft import fft_convolve as sp_fft_convolve

P_MIN = 1e-3
KERNEL_CUTOFF = 600

def convolve(signal, sigma, pMin=None, cutoff=None):
	if cutoff is None: cutoff = KERNEL_CUTOFF

	kernel = convolution_kernel(sigma, len(signal), pMin)	

	if len(kernel) < cutoff:
		convSignal = convolve1d(signal, kernel, mode='constant', cval=0.0)
	else:
		convSignal = fftw_convolve(signal, kernel)
	return convSignal
	
def convolution_kernel(sigma, signalLength, pMin=None):
	if pMin is None: pMin = P_MIN

	normDistr = norm(loc=0, scale=sigma)
	x = round(normDistr.ppf(pMin))
	kernel = normDistr.pdf(arange(x, (-1 * x) + 1))

	return kernel / kernel.sum()