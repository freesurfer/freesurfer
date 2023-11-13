#!/usr/bin/env python

import re
import pathlib
import setuptools

# extract the current version
init_file = pathlib.Path(__file__).parent.resolve().joinpath('spheremorph/__init__.py')
init_text = open(init_file, 'rt').read()
pattern = r"^__version__ = ['\"]([^'\"]*)['\"]"
match = re.search(pattern, init_text, re.M)
if not match:
    raise RuntimeError(f'Unable to find __version__ in {init_file}')
version = match.group(1)

# run setup
setuptools.setup(
    name='spheremorph',
    version=version,
    license='GPLv3',
    description='Spherical Registration with Convolutional Networks',
    url='https://github.com/silencer1127/spheremorph',
    keywords=['deformation', 'registration', 'imaging', 'cnn', 'mri'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'License :: GNU General Public License v3',
        'Operating System :: OS Independent',
    ],
    install_requires=[
        'numpy',
        'tensorflow>=2.0.0',
        'keras',
        'voxelmorph',
        'neurite',
    ]
)
