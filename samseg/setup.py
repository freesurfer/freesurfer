#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os

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
    'pyqtgraph',
    'pyqt5',
]

here = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = '\n' + f.read()

our_packages = find_packages(include=['samseg', 'samseg.dev_utils', 'samseg.hdav.hdav'])

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
    # $ setup.py publish support.
    # cmdclass={
    #     'upload': UploadCommand,
    # },
)
