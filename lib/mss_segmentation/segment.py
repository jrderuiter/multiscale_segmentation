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

from math import log, ceil
from numpy import asarray, nonzero, power, exp, sqrt, arange, diff, unique, append, mean
from numpy import int as int_type
from numpy import double as np_double

from convolution import convolve
from linking.link import link, LINK_TYPES

DELTA_TAU = 0.5 * log(2)

def segment_signal(signal, numScales, deltaTau=None, kMin=None, linkType=None):
	if deltaTau is None: deltaTau = DELTA_TAU
	if kMin is None: kMin = pow((1 - exp(-2 * deltaTau)), -0.5)

	signal = asarray(signal).flatten()
	signalLength, maxDiff, numLeadingZeros, numTrailingZeros = _signal_properties(signal)

	# Calculate scale sigma's
	scaleSigmaArray = exp(arange(0, numScales) * deltaTau)

	# Setup initial node set
	nodeIds = nonzero(diff(signal, axis=0))[0] + 1
	if nodeIds[-1] != signalLength: 
		nodeIds = append(nodeIds, signalLength)
	
	# Loop!
	baseSignal = asarray(signal, dtype=np_double)
	prevScaleSignal, prevSegmentEnds = baseSignal, nodeIds
	nodeMapping, segmentEnds  = {}, { 0: nodeIds }
	for scaleIndex in range(1, numScales):
		scaleSignal = convolve(baseSignal, scaleSigmaArray[scaleIndex])

		d, dcp, r = _search_volume(scaleSigmaArray, scaleIndex, kMin)
		
		parents, scaleSegmentEnds = link(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, \
										 maxDiff, prevSegmentEnds, numLeadingZeros, numTrailingZeros, linkType=linkType)	

		nodeMapping[scaleIndex] = zip(nodeIds, parents)
		segmentEnds[scaleIndex] = scaleSegmentEnds

		prevScaleSignal, prevSegmentEnds, nodeIds = scaleSignal, scaleSegmentEnds, unique(parents)

	return nodeMapping, segmentEnds

def segment_ends_to_segments(segmentEnds):
	segments = {}
	for i, endArray in segmentEnds.items():
		segments[i] = zip(endArray[:-1], endArray[1:])
		segments[i].insert(0, (0, endArray[0]))
	return segments

def _signal_properties(signal):
	signalLength = len(signal)
	maxDiff = signal.max() - signal.min()
	
	nonZeroIndices = nonzero(signal)[0]
	numLeadingZeros = nonZeroIndices[0]
	numTrailingZeros = signalLength - nonZeroIndices[-1] - 1
	
	return signalLength, maxDiff, numLeadingZeros, numTrailingZeros

def _search_volume(scaleSigmaArray, scaleIndex, kMin):
	currSigma, prevSigma = scaleSigmaArray[scaleIndex], scaleSigmaArray[scaleIndex-1]

	r = scaleSigmaArray[scaleIndex] if scaleIndex == 1 else sqrt(pow(currSigma, 2) - pow(prevSigma, 2)) 
	r = ceil(kMin * r)
		
	dcp = arange(-1 * r, r + 1, dtype=int_type)
	
	d = exp(power(dcp, 2) / (-2 * (pow(currSigma, 2) - pow(prevSigma, 2))))
	d /= exp(pow(0.5 * currSigma, 2) / (-2 * (pow(currSigma, 2) - pow(prevSigma, 2))))
			
	j = nonzero(abs(dcp) <= 0.5 * currSigma)[0]
	d[j] = 1 + 1e-6 * d[j]
	
	return d, dcp, r		# Transpose d for later calculations
