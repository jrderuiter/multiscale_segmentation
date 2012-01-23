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

from numpy import unique, subtract, clip, multiply, abs, argmax, empty
from numpy import int as int_type
from .py_link_shared import calc_parents_and_segment_ends, calc_gp, W1

import pdb

def link_parents(nodeIds, d, dcp, r, prevScaleSignal, scaleSignal, signalLength, maxDiff, prevSegmentEnds, nLZ, nTZ):
	uk, px1, px2, x1, x2 = _init_pass(signalLength, dcp, r, nodeIds)
	parents, scores = empty(len(uk), dtype=int_type), []
	d = d.flatten()

	# Pass 1
	for i in range(len(uk)):
		score = multiply(d[x1[i]:x2[i]], (1 - (abs(prevScaleSignal[uk[i]] - scaleSignal[px1[i]:px2[i]])) / maxDiff))
		parents[i] = argmax(score)
		scores.append(score)
	parents += px1	
	parents, scaleSegmentEnds = calc_parents_and_segment_ends(parents, prevSegmentEnds, signalLength)

	# Calc GP
	gp = calc_gp(signalLength, scaleSegmentEnds, parents, nLZ, nTZ)

	# Pass 2
	for i in range(len(uk)):
		scores[i] += d[x1[i]:x2[i]] * ((W1 * gp[px1[i]:px2[i]]) / gp.max())
		parents[i] = argmax(scores[i])
	parents += px1	
	parents, scaleSegmentEnds = calc_parents_and_segment_ends(parents, prevSegmentEnds, signalLength)

	return parents, scaleSegmentEnds

def _init_pass(signalLength, dcp, r, nodeIds):
	uk = unique(nodeIds)
	
	px1, px2 = uk + dcp[0], uk + dcp[-1] + 1
	clip(px1, 0, signalLength, px1)
	clip(px2, 0, signalLength, px2)
	
	x1 = uk - px1
	subtract(r, x1, x1)
	
	x2 = px2 - uk
	subtract(r, x2, x2)
	subtract(len(dcp)-1, x2, x2)

	return uk, px1, px2, x1, x2