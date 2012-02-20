from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

cFiles =  ['lib/mss_segmentation/linking/c_ext/link.c', 'lib/mss_segmentation/linking/c_ext/link_window.c', 
		   'lib/mss_segmentation/linking/c_ext/link_parents.c', 'lib/mss_segmentation/linking/c_ext/link_shared.c']

setup(
	name = 'mss-segmentation',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('lib.mss_segmentation.linking.c_ext.link', cFiles, extra_compile_args=['-std=c99'])]
)
