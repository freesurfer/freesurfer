#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
from setuptools import setup, find_packages, Distribution


class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

    def has_ext_modules(self):
        return True


def find_gems_libraries():
    search_pattern = os.path.join('gems', 'gems.*.so')
    libraries = [os.path.basename(filename) for filename in glob.glob(search_pattern)]
    if not libraries:
        print('error: could not find gems library')
        sys.exit(1)
    return libraries


setup(
    distclass=BinaryDistribution,
    name='gems',
    version='0.1.0',
    description="Python wrapper for GEMS2 C++ tooling",
    package_dir={'gems': 'gems'},
    packages=['gems'],
    package_data={'gems': find_gems_libraries()},
    zip_safe=False,
    keywords='samseg gems freesurfer',
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
    ],
)
