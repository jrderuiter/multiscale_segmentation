from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('mss_segmentation.linking.c_ext.link', ['mss_segmentation/linking/c_ext/link.c', 'mss_segmentation/linking/c_ext/link_window.c', 'mss_segmentation/linking/c_ext/link_parents.c', 'mss_segmentation/linking/c_ext/link_shared.c'])]
)