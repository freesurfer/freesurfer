#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys

from setuptools import setup, find_packages

# Package meta-data.
NAME = 'samseg'
DESCRIPTION = 'GEMS2 Samseg.'
URL = 'https://github.com/innolitics/freesurfer'
EMAIL = 'wmills@innolitics.com, yshrestha@innolitics.com'
AUTHOR = 'William Mills, Yujan Shrestha'
VERSION = '0.1.1'

REQUIRED = [
    'nibabel >= 2.2.1',
    'numpy >= 1.13.3',
    'scipy >= 1.0.0',
    'gems2python',
]

PACKAGES_TO_FIND = ['samseg', 'samseg.dev_utils']
NO_VISUALIZATION_ARG = '--no-visualization'
if not NO_VISUALIZATION_ARG in sys.argv:
    REQUIRED += [ 'pyqtgraph', 'pyqt5']
    PACKAGES_TO_FIND += ['samseg.hdav.hdav']
else:
    sys.argv.remove(NO_VISUALIZATION_ARG)

here = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = '\n' + f.read()

our_packages = find_packages(include=PACKAGES_TO_FIND)

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    packages=our_packages,
    install_requires=REQUIRED,
    include_package_data=True,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    entry_points={
        'console_scripts': [
            'run_samseg = samseg.main:samseg_main',
        ]
    }
)
