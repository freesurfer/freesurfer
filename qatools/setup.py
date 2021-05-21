from setuptools import setup, find_packages
from os import path

with open("README.md", "r") as fh:
    long_description = fh.read()

# see: https://packaging.python.org/guides/distributing-packages-using-setuptools/#setup-py for a description
# see: https://github.com/pypa/sampleproject/blob/master/setup.py for an example
# see: https://pip.pypa.io/en/stable/reference/pip_install for a description how to use pip install in conjunction with a git repository

setup(
    name="qatoolspython",
    version="0.9.6-beta",
    author="Kersten Diers, Martin Reuter, and others (see README)",
    description="A set of quality assurance / quality control scripts for Freesurfer 6.0 processed structural MRI data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/reuter-lab/qatools-python",
    packages=['qatoolspython'],
    # see https://pypi.org/classifiers/
    classifiers=[ 
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT",
#        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
    keywords='Freesurfer',
)
