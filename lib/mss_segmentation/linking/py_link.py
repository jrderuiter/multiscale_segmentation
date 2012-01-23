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

from .py_link_parents import link_parents
from .py_link_window import link_window

CUT_OFF = 300

def link_composite(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, segmentEnds, nLZ, nTZ):
	if len(nodeIds) > CUT_OFF:
		result = link_window(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, segmentEnds, nLZ, nTZ)
	else:
		result = link_parents(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, segmentEnds, nLZ, nTZ)
	return result