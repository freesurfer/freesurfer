#!/usr/bin/env python

import sys
import glob
import platform
import operator
import os.path as path
from setuptools import setup, find_packages, Distribution


# the freesurfer python packages
packages = [
    'freesurfer',
    'freesurfer.utils',
    'freesurfer.samseg',
    'freesurfer.deeplearn',
    'freesurfer.subregions',
]

# get required dependencies from requirements.txt
with open(path.join(path.dirname(path.realpath(__file__)), 'requirements.txt')) as file:
    requirements = [line for line in file.read().splitlines() if not line.startswith('#')]

# ---- run the setup ----

# since the freesurfer package will contain python-wrapped c++ libraries, we need to indicate
# that it will be a platform-dependent distribution via a custom Distribution class
class BinaryDistribution(Distribution):
    def is_pure(self):
        return False
    def has_ext_modules(self):
        return True

# locates cpython libraries compiled with pybind
def find_libs(libname, required=True):
    libraries = glob.glob('**/%s.*%s*.so' % (libname, platform.system().lower()), recursive=True)
    if required and not libraries:
        print('error: could not find %s library that matches the current python version' % libname)
        sys.exit(1)
    return [path.basename(filename) for filename in libraries]

setup(
    distclass=BinaryDistribution,
    name='freesurfer',
    version='0.0.1',
    description='Python package for FreeSurfer neuroimaging software',
    author='Laboratory for Computational Neuroimaging',
    author_email='freesurfer@nmr.mgh.harvard.edu',
    url='https://github.com/freesurfer/freesurfer',
    packages=find_packages(include=packages),
    package_data={
        'freesurfer': operator.add(
            find_libs('bindings'),
            find_libs('labelfusion', required=False)),
        'freesurfer.samseg':
            find_libs('gemsbindings')
    },
    install_requires=requirements,
    include_package_data=True
)
