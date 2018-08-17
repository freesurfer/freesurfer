#!/usr/bin/env python

import os
import sys
import glob
from setuptools import setup, find_packages, Distribution
from distutils.util import convert_path


requirements = [
    'nibabel == 2.2.1',
    'numpy >= 1.13.3',
    'scipy >= 1.0.0'
]

packages = [
    'freesurfer',
    'freesurfer.gems'
]

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False
    def has_ext_modules(self):
        return True

def find_shared_libraries(libname):
    libraries = glob.glob('**/%s.*.so' % libname, recursive=True)
    if not libraries:
        print('error: could not find %s library' % libname)
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
    package_data={'freesurfer.gems': find_shared_libraries('gems')},
    install_requires=requirements,
    include_package_data=True
)
