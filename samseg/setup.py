#!/usr/bin/env python

import io
import os
import sys
from setuptools import setup, find_packages


REQUIRED = [
    'nibabel == 2.2.1',
    'numpy >= 1.13.3',
    'scipy >= 1.0.0'
]

PACKAGES_TO_FIND = ['samseg', 'samseg.dev_utils']
NO_VISUALIZATION_ARG = '--no-visualization'
if NO_VISUALIZATION_ARG not in sys.argv:
    REQUIRED += ['pyqtgraph', 'pyqt5']
    PACKAGES_TO_FIND += ['samseg.hdav.hdav']
else:
    sys.argv.remove(NO_VISUALIZATION_ARG)

setup(
    name='samseg',
    version='0.1.1',
    description='sequence-adaptive multimodal segmentation',
    author='William Mills, Yujan Shrestha',
    author_email='wmills@innolitics.com, yshrestha@innolitics.com',
    url='https://github.com/innolitics/freesurfer',
    packages=find_packages(include=PACKAGES_TO_FIND),
    install_requires=REQUIRED,
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: Implementation :: CPython',
    ]
)
