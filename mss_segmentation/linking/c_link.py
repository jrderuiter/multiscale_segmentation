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

from numpy.ctypeslib import ndpointer, load_library
from numpy import asarray, empty, unique
from numpy import intc as np_intc
from numpy import double as np_double
from ctypes import c_int, c_double
from os import path

CUT_OFF = 300

_link = load_library('liblink', path.join(path.dirname(__file__), 'c_ext'))

_link.link_parents.argtypes = [ndpointer(dtype=np_intc), c_int, c_double, ndpointer(dtype=np_double), ndpointer(dtype=np_intc), c_int,\
							   ndpointer(dtype=np_double), ndpointer(dtype=np_double), c_int, c_double, ndpointer(dtype=np_intc), c_int,
					   		   c_int, c_int, ndpointer(dtype=np_intc), ndpointer(dtype=np_intc)]
_link.link_parents.restype = c_int

_link.link_window.argtypes = [ndpointer(dtype=np_intc), c_int, ndpointer(dtype=np_double), ndpointer(dtype=np_intc), c_int,\
							   ndpointer(dtype=np_double), ndpointer(dtype=np_double), c_int, c_double, ndpointer(dtype=np_intc), c_int,
					   		   c_int, c_int, ndpointer(dtype=np_intc), ndpointer(dtype=np_intc)]
_link.link_window.restype = c_int

def link_composite(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ):
	if len(nodeIds) > CUT_OFF:
		result = link_window(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ)
	else:
		result = link_parents(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ)
	return result

def link_parents(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ):
	cUk, cDcp, cPrevSegmentEnds = _setup_c_vars(nodeIds, dcp, prevSegmentEnds)

	lenUk, lenSegments = len(cUk), len(prevSegmentEnds)
	parents, segmentEnds = empty(lenUk, dtype=np_intc), empty(lenSegments, dtype=np_intc)

	numSegments = _link.link_parents(cUk, lenUk, r, d, cDcp, len(dcp), prevScaleSignal, scaleSignal, signalLength, maxDiff, \
							 		 cPrevSegmentEnds, lenSegments, nLZ, nTZ, parents, segmentEnds);

	return parents, segmentEnds[:numSegments]

def link_window(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ):
	cUk, cDcp, cPrevSegmentEnds = _setup_c_vars(nodeIds, dcp, prevSegmentEnds)

	lenUk, lenSegments = len(cUk), len(prevSegmentEnds)
	parents, segmentEnds = empty(lenUk, dtype=np_intc), empty(lenSegments, dtype=np_intc)

	numSegments = _link.link_window(cUk, lenUk, d, cDcp, len(dcp), prevScaleSignal, scaleSignal, signalLength, maxDiff, \
							 		cPrevSegmentEnds, lenSegments, nLZ, nTZ, parents, segmentEnds);

	return parents, segmentEnds[:numSegments]

def _setup_c_vars(nodeIds, dcp, prevSegmentEnds):
	cUk = asarray(unique(nodeIds), dtype=np_intc)
	cDcp = asarray(dcp, dtype=np_intc)
	cPrevSegmentEnds = asarray(prevSegmentEnds, dtype=np_intc)

	return cUk, cDcp, cPrevSegmentEnds
