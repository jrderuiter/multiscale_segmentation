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

from numpy import nonzero, diff, append, zeros
import pdb

W1 = 1

def calc_parents_and_segment_ends(parents, prevSegmentEnds, signalLength):
	parents.sort()

	uniqueParents = nonzero(diff(parents))
	scaleSegmentEnds = prevSegmentEnds[uniqueParents]
	scaleSegmentEnds = append(scaleSegmentEnds, signalLength)
	
	return parents, scaleSegmentEnds

def calc_gp(signalLength, scaleSegmentEnds, parents, numLeadingZeros, numTrailingZeros):
	uniqueParents = append(parents[nonzero(diff(parents))], parents[-1])
		
	gp = zeros(signalLength)
	gp[uniqueParents] = append(scaleSegmentEnds[0] + 1, diff(scaleSegmentEnds))
	
	#gp[uniqueParents[0]] -= numLeadingZeros     # Check if necessary in algorithm?
	#gp[uniqueParents[-1]] -= numTrailingZeros	# Seems to affect conversion, disabling for now.
	
	#pdb.set_trace()

	return gp
