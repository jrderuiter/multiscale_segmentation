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

from numpy import unique, repeat, logical_and, abs, argmax, newaxis
from .py_link_shared import calc_gp, calc_parents_and_segment_ends, W1

def link_window(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ):
	uk = unique(nodeIds)
	lenDcp = len(dcp)

	ukRepeat = repeat(uk[newaxis,:], lenDcp, 0)
	ukDcp = uk + dcp.reshape((lenDcp, 1))
	ukMask = logical_and(ukDcp >= 0, ukDcp < signalLength)
	ukDcp *= ukMask

	d = d.reshape((len(d),1))

	# Pass 1
	score = d * (1 - (abs(prevScaleSignal[ukRepeat] - scaleSignal[ukDcp]) / maxDiff))
	parents, scaleSegmentEnds = _conclude_pass(uk, ukMask, signalLength, dcp, prevSegmentEnds, score)

	# Calc GP
	gp = calc_gp(signalLength, scaleSegmentEnds, parents, nLZ, nTZ)

	# Pass 2
	score += d * ((W1 * gp[ukDcp]) / gp.max())
	parents, scaleSegmentEnds = _conclude_pass(uk, ukMask, signalLength, dcp, prevSegmentEnds, score)

	return parents, scaleSegmentEnds
	
def _conclude_pass(uk, ukMask, signalLength, dcp, prevSegmentEnds, score):
	score *= ukMask
	parents = uk + dcp[argmax(score, axis=0)]
	
	return calc_parents_and_segment_ends(parents, prevSegmentEnds, signalLength)
	