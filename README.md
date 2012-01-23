Multiscale Segmentation: Scale-space segmentation of one-dimensional signals
=============================================================================

## Description

This is a python library that implements a scale-space segmentation algorithm very similar to that described in the following paper: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=574787&tag=1. The implementation is currently limited to one-dimensional signals, though the code should be easily extendable.

The algorithm employs three main steps in its segmentation procedure:

1. Blurring the signal to create a scale-space signal representation.
2. Building a segmentation tree across the scale-space using a linking algorithm.
3. Segmenting the scale-space using the segmentation tree.

Step 1 is implemented using convolve1d and fftconvolve from the scipy software package, where convolve1d is used for small inputs and fftconvolve for larger inputs to speed up the processing. Step 2 is implemented in both python and C (for more speed). The C implementation is used by default, though this does require building of the C extension files.

## Installation

The library can currently be used by adding the libary path to your python path. Note that the C extentions must be built first, which can be done with the setup.py script using the following command in the library directory: 

python setup.py build_ext --inplace

## Usage

Currently, a signal can be segmented by passing a vector of its values to the segment_signal function along with the number of desired scales. This returns a mapping dictionary representing the build segmentation tree and a dictionary of segment ends, which contains the end of each segment per scale. The dictionary of segment ends can be transformed to segment tuples (start, end) using the segment_ends_to_segments function.

The expansion module provides two convenience methods for imposing the segmentation on the original signal. By passing the segment ends dictionary to the expand_segmentation function along with the original signal, a two-dimensional matrix is calculated that represents the segmented multi-scale representation of the signal.

## Examples

Currently a single example script (signal_to_multiscale.py) is included in the examples directory. This script allows the user to specify an input signal (currently in the form of a matlab file containing a signal), an output file (also a matlab file) and the number of scales with which the multi-scale analysis is to be performed. The script reads the input signal, applies the multi-scale segmentation algorithm and saves the output as a 2D matrix in the specified output file (which is created). This file can be opened for futher analysis, for example by loading it in matlab and visualizing the result with image or imagesc. Additional file types are easy to add, as long as the input signal can be converted to a 1D list or array and the output file can store 2D matrices.

Usage: python ./scripts/signal_to_multiscale.py --input ./data/testSignal.mat --output output.mat --scales 35