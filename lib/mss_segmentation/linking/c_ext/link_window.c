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

#include "link_window.h"

#define W1 1

extern void calc_gp(int*, int, int*, int, int, int, int, double*);
extern double max(double*, int);
extern int compare_int(const void*, const void*);

static void link_pass1(int*, int, double*, int*, int, double*, double*, int, double, int *);
static void link_pass2(int*, int, double*, int*, int, double*, double*, int, double, double*, int, int*);
static int conclude_pass(int*, int, int, int*, int, int*);

int link_window(int *uk, int lenUK, double *d, int *dcp, int lenDCP, \
		 double *prevScaleSignal, double *scaleSignal, int signalLength, double diffMax, \
		 int* prevSegmentEnds, int lenPrevSegEnds, int nTZ, int nEZ, int *parents, int *segmentEnds) {

	int lenSegmentEnds;
	double *gp = malloc(sizeof(double) * signalLength);

	link_pass1(uk, lenUK, d, dcp, lenDCP, prevScaleSignal, scaleSignal, signalLength, diffMax, parents);
	lenSegmentEnds = conclude_pass(parents, lenUK, signalLength, prevSegmentEnds, lenPrevSegEnds, segmentEnds);
	
	calc_gp(parents, lenUK, segmentEnds, lenSegmentEnds, signalLength, nTZ, nEZ, gp);

	link_pass2(uk, lenUK, d, dcp, lenDCP, prevScaleSignal, scaleSignal, signalLength, diffMax, gp, signalLength, parents);
	lenSegmentEnds = conclude_pass(parents, lenUK, signalLength, prevSegmentEnds, lenPrevSegEnds, segmentEnds);

	free(gp);
	return lenSegmentEnds;
}

static void link_pass1(int *uk, int lenUK, double *d, int *dcp, int lenDCP, \
				double *prevScaleSignal, double *scaleSignal, int lenSignal, double diffMax, int *parents) {
	int newPos, i, j, bestScoringIndex;
	double score, bestScore;

	for(j=0; j < lenUK; j++) {
		score = 0;
		bestScore = 0;

		for(i=0; i < lenDCP; i++) {
			newPos = uk[j] + dcp[i];
			if(newPos >= 0 && newPos < lenSignal) {	
				score = d[i] * (1 - (fabs(prevScaleSignal[uk[j]] - scaleSignal[newPos])/diffMax));
				if(score > bestScore) {
					bestScore = score;
					bestScoringIndex = i;
				}
			}
		}

		parents[j] = uk[j] + dcp[bestScoringIndex];
	}
}

static void link_pass2(int *uk, int lenUK, double *d, int *dcp, int lenDCP, \
				double *prevScaleSignal, double *scaleSignal, int lenSignal, double diffMax, \
				double *gp, int lenGP, int *parents) {
	int newPos, i, j, bestScoringIndex;
	double score, bestScore, 
		   gpMax = max(gp, lenGP);

	for(j=0; j < lenUK; j++) {
		score = 0;
		bestScore = 0;

		for(i=0; i < lenDCP; i++) {
			newPos = uk[j] + dcp[i];
			if(newPos >= 0 && newPos < lenSignal) {
				score = d[i] * ((1 - (fabs(prevScaleSignal[uk[j]] - scaleSignal[newPos])/diffMax)) + \
							    (W1 * gp[newPos]/gpMax));
				if(score > bestScore) {
					bestScore = score;
					bestScoringIndex = i;
				}
			}
		}

		parents[j] = uk[j] + dcp[bestScoringIndex];
	}
}

static int conclude_pass(int *parents, int lenParents, int signalLength, int *prevSegmentEnds, int lenPrevSegEnds, int *segmentEnds) {	
	qsort(parents, lenParents, sizeof(int), compare_int);

	int i, segmentCounter = 0;
	for(i=0; i < lenParents-1; i++)
		if(parents[i] < parents[i+1]) 
			segmentEnds[segmentCounter++] = prevSegmentEnds[i];
	segmentEnds[segmentCounter++] = signalLength;

	return segmentCounter;	// Length segments
}

