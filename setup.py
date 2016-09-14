#!/usr/bin/env python
import os
from distutils.core import setup
from distutils.extension import Extension
import glob
import re

import numpy

if os.path.exists("bsplines.pyx"):
    USE_CYTHON = True
    fname = "bsplines.pyx"
else:
    USE_CYTHON = False
    fname = "bsplines.c"

sourcefiles = [fname, os.path.join("src", "bspl.c")]
headerfiles = [os.path.join("src", "bspl.h")]
include_dirs=[numpy.get_include(), "src"]
extensions = [Extension("bsplines", sourcefiles, include_dirs=include_dirs,
                        depends=headerfiles)]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

# Synchronize version from code.
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

classifiers = []

setup(name="bsplines", 
      version=version,
      description="",
      long_description="",
      license="MIT",
      classifiers=classifiers,
      url="",
      author="Kyle Barbary",
      author_email="kylebarbary@gmail.com",
      ext_modules=extensions)
