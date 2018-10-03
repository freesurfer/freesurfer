#!/usr/bin/env python

import os
import sys
import glob
import platform
from setuptools import setup, find_packages, Distribution


# the freesurfer python packages
packages = [
    'freesurfer',
    'freesurfer.gems',
    'freesurfer.samseg'
]

# required dependencies
requirements = [
    'nibabel == 2.3.0',
    'numpy == 1.13.3',
    'scipy == 1.0.0'
]


# ---- run the setup ----

# since the freesurfer package will contain python-wrapped c++ libraries, we need to indicate
# that it will be a platform-dependent distribution via a custom Distribution class
class BinaryDistribution(Distribution):
    def is_pure(self):
        return False
    def has_ext_modules(self):
        return True

# locates cpython libraries compiled with pybind
def find_shared_libs(libname):
    libraries = glob.glob('**/%s.*%s*.so' % (libname, platform.system().lower()), recursive=True)
    if not libraries:
        print('error: could not find %s library that matches the current python version' % libname)
        sys.exit(1)
    return [os.path.basename(filename) for filename in libraries]

setup(
    distclass=BinaryDistribution,
    name='freesurfer',
    version='0.0.1',
    description='Python package for FreeSurfer neuroimaging software',
    author='Laboratory for Computational Neuroimaging',
    author_email='freesurfer@nmr.mgh.harvard.edu',
    url='https://github.com/freesurfer/freesurfer',
    packages=find_packages(include=packages),
    package_data={'freesurfer.gems': find_shared_libs('gems_python')},
    install_requires=requirements,
    include_package_data=True
)
