from pyfft.cl import Plan
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from pyopencl.elementwise import ElementwiseKernel

def fft_convolve(signal, kernel):
	ctx = cl.create_some_context(interactive=True)
	queue = cl.CommandQueue(ctx)

	s1, s2 = len(signal), len(kernel)
	size = s1 + s2 - 1
	fsize = int(2 ** np.ceil(np.log2(size)))

	paddedSig = _zero_pad_array(signal, fsize, dtype=np.complex64)
	paddedKern = _zero_pad_array(kernel, fsize, dtype=np.complex64)

	plan = Plan((fsize,1), queue=queue, dtype=np.complex64, fast_math=False)
	gpu_sig = cl_array.to_device(ctx, queue, paddedSig)
	gpu_kern = cl_array.to_device(ctx, queue, paddedKern)

	plan.execute(gpu_sig.data)
	plan.execute(gpu_kern.data)

	multKern = ElementwiseKernel(ctx, "cfloat_t *sig, cfloat_t *kern", "sig[i] = cfloat_mul(sig[i], kern[i])")
	multKern(gpu_sig, gpu_kern)

	plan.execute(gpu_sig.data, inverse=True)
	result = gpu_sig.get()

	fslice = slice((s2-1)/2, (s2-1)/2+s1)
	return np.real(result[fslice])