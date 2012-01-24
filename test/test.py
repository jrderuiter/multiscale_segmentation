import sys
from os import path
sys.path.insert(0, path.abspath(path.join(path.dirname(__file__), '../lib')))

from mss_segmentation.segment import segment_signal, LINK_TYPES

from pprint import pprint as pp
from scipy.io import loadmat

def main(argv):
	refSignal = loadmat('refSignal')['V'][:,0]
	refMapping = refKM_to_mapping(loadmat('refKM')['KM'])

	segmentEnds, mapping = segment_signal(refSignal, len(refMapping), linkType=LINK_TYPES.C_COMPOSITE, doNodeMapping=True)

	print 'Comparing to reference'
	for i in range(1, len(refMapping)):
		print "- Scale %d equal: %s" % (i, str(refMapping[i] == mapping[i]))

def refKM_to_mapping(refKM):
	mapping = {}

	m, n = refKM.shape
	for i in range(1, m):
		kids = refKM[i][0][0] - 1 
		parents = refKM[i][0][1] - 1

		mapping[i] = zip(kids, parents)

	return mapping

def difference(refMap, currMap):
	return [tup for tup in refMap if tup not in currMap]

if __name__ == '__main__':
	main(sys.argv)