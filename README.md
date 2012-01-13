Multiscale Segmentation - Scale-space segmentation of one-dimensional signals
=============================================================================

## Description

This is a python library that implements a scale-space segmentation algorithm very similar to that described in the following paper: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=574787&tag=1. The implementation is currently limited to one-dimensional signals, though the code should be easily extendable.

The algorithm employs three main steps in its segmentation procedure:

1. Blurring the signal to create a scale-space signal representation.
2. Building a segmentation tree across the scale-space using a linking algorithm.
3. Segmenting the scale-space using the segmentation tree.

Step 1 is implemented using convolve1d and fftconvolve from the scipy software package, where convolve1d is used for small inputs and fftconvolve for larger inputs to speed up the processing. Step 2 is implemented in both python and C (for more speed). The C implementation is used by default, though this does require building of the C extension files.

## Usage

Currently, a signal can be segmented by passing a vector of its values to the segment_signal function along with the number of desired scales. This returns a mapping dictionary representing the build segmentation tree and a dictionary of segment ends, which contains the end of each segment per scale. The dictionary of segment ends can be transformed to segment tuples (start, end) using the segment_ends_to_segments function.

The expansion module provides two convenience methods for imposing the segmentation on the original signal. By passing the segment ends dictionary to the expand_segmentation function along with the original signal, a two-dimensional matrix is calculated that represents the segmented multi-scale representation of the signal.

Note that this package will be expanded with a number of (example) scripts that allow a user to specify input- and output-files as a command-line option, thus not requiring any programming for use of this library.