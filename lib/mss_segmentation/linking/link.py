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

import c_link
import py_link

def enum(*sequential, **named):
	enums = dict(zip(sequential, range(len(sequential))), **named)
	return type('Enum', (), enums)

LINK_TYPES = enum('C_PARENTS', 'C_WINDOW', 'C_COMPOSITE', 'PY_PARENTS', 'PY_WINDOW', 'PY_COMPOSITE')

funcDict = {
	LINK_TYPES.C_PARENTS: c_link.link_parents,
	LINK_TYPES.C_WINDOW: c_link.link_window,
	LINK_TYPES.C_COMPOSITE: c_link.link_composite,
	LINK_TYPES.PY_PARENTS: py_link.link_parents,
	LINK_TYPES.PY_WINDOW: py_link.link_window,
	LINK_TYPES.PY_COMPOSITE: py_link.link_composite
}

def link(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ, linkType=None):
	if linkType is None: linkType = LINK_TYPES.C_COMPOSITE
	assert linkType in funcDict
	
	return funcDict[linkType](nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ)