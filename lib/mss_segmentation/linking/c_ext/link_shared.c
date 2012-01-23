/*
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
*/

#include <stdio.h>

double max(double *array, int lenArray) {
	double maxValue = -9999999;

	for(int i=0; i < lenArray; i++)
		if(maxValue < array[i])
			maxValue = array[i];
	
	return maxValue;
}

int compare_int(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

/*
UP = [P(diff(P)~=0) P(end)];
Gp = zeros(N,1);
Gp(UP) = [SegmentEnd{n}(1) diff(SegmentEnd{n})];
Gp(UP(1)) = Gp(UP(1))-nTZ;
Gp(UP(end)) = Gp(UP(end))-nEZ;
*/

void calc_gp(int *parents, int parentLen, int* segmentEnds, int segmentEndsLength, int signalLength, int nTZ, int nEZ, double *gp) {
	// Fill in GP
	int i, nextPIndex,
		upCounter = -1, pIndex = 0;

	while(pIndex < parentLen-1 && parents[pIndex] == parents[pIndex+1])	 // Stop if pIndex is unique or last parent
		pIndex++;

	for(i=0; i < signalLength; i++) {
		if(pIndex < parentLen && i == parents[pIndex]) {
			nextPIndex = pIndex + 1;
			while(nextPIndex < parentLen-1 && parents[nextPIndex] == parents[nextPIndex+1])
				nextPIndex++;

			if(upCounter == -1)	// -1 signals first item
				gp[i] = segmentEnds[0] + 1;
			else
				gp[i] = segmentEnds[upCounter+1] - segmentEnds[upCounter];

			pIndex = nextPIndex;
			upCounter++;
		} else
			gp[i] = 0;
	}

	//gp[parents[parentLen-1]] -= nEZ; 
}
