
import sys
from os import path
import argparse

from scipy.io import loadmat, savemat

sys.path.insert(0, path.abspath(path.join(path.dirname(__file__), '../mss_segmentation')))

from segment import segment_signal
from expansion import expand_segmentation

def main(argv):
	args = parse_args(argv)

	inputSignal = open_input(args.input, args.signal)
	
	mapping, segmentEnds = segment_signal(inputSignal, args.scales)
	expandedSegmentation = expand_segmentation(inputSignal, segmentEnds)

	save_output(args.output, expandedSegmentation)

def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help='Input file path', required=True)
	parser.add_argument('--signal', help='Name of signal in input file', default='signal')
	parser.add_argument('--output', help='Output file path', required=True)
	parser.add_argument('--scales', help='Number of scales', required=True, type=int)

	return parser.parse_args()

def open_input(filePath, signalName):
	if path.splitext(filePath)[1] == '.mat':
		return open_matlab_input(filePath, signalName)
	else:
		raise Exception('Unsupported output file type, currently only supports Matlab (*.mat) input files.')

def open_matlab_input(filePath, signal_name):
	return loadmat(filePath)[signal_name]

def save_output(filePath, outputSignal):
	if path.splitext(filePath)[1] == '.mat':
		return save_matlab_output(filePath, outputSignal)
	else:
		raise Exception('Unsupported output file type, currently only supports Matlab (*.mat) ouput files.')

def save_matlab_output(filePath, outputSignal):
	savemat(filePath, { 'signal': outputSignal }, do_compression=True, oned_as='row')

if __name__ == '__main__':
	main(sys.argv)