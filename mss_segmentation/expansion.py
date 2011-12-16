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

from numpy import vstack, empty, mean, array, searchsorted

def expand_segmentation(originalSignal, segmentEnds):
	if type(segmentEnds) is dict:
		expandedSignal = vstack([expand_segmentation(originalSignal, segmentEnds[scale]) for scale in range(len(segmentEnds))])
	else:
		expandedSignal = empty(len(originalSignal))

		segmentStart = 0
		for segmentEnd in segmentEnds:
			expandedSignal[segmentStart:segmentEnd] = mean(originalSignal[segmentStart:segmentEnd])
			segmentStart = segmentEnd

	return expandedSignal

def expand_segmentation_for_range(originalSignal, segmentEnds, rangeTuple):
	assert rangeTuple[0] >= 0 and rangeTuple[0] < len(originalSignal)
	assert rangeTuple[1] > 0 and rangeTuple[1] <= len(originalSignal)

	if type(segmentEnds) is dict:
		expandedSignal = vstack([expand_segmentation_for_range(originalSignal, segmentEnds[scale], rangeTuple) for scale in range(len(segmentEnds))])
	else:
		rangeStart, rangeEnd = rangeTuple
		firstSegIndex = searchsorted(segmentEnds, rangeStart)
		lastSegIndex = searchsorted(segmentEnds, rangeEnd)

		expandedSignal = empty(rangeEnd - rangeStart)

		startIndex = 0
		segmentStart = segmentEnds[firstSegIndex-1] if firstSegIndex > 0 else 0
		for segmentEnd in segmentEnds[firstSegIndex:lastSegIndex+1]:
			endIndex = min(segmentEnd - rangeStart, rangeEnd)
			expandedSignal[startIndex:endIndex] = mean(originalSignal[segmentStart:segmentEnd])
			segmentStart, startIndex = segmentEnd, endIndex
		
	return expandedSignal
