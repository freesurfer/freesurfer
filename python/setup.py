#!/usr/bin/env python

import sys
import glob
import platform
import operator
import os
import os.path as path
from setuptools import setup, find_packages, Distribution

# the freesurfer python packages
INTEGRATE_SAMSEG = os.environ.get("INTEGRATE_SAMSEG_OPTION")
if (INTEGRATE_SAMSEG == 'OFF'):
  packages = [
    'fsbindings',
    'gems',
    'gems.subregions',
  ]
else:
  packages = [
    'fsbindings',
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
    libraries = glob.glob('%s.*%s*.so' % (libname, platform.system().lower()), recursive=True)
    if not libraries:
        libraries = glob.glob('**/%s.*%s*.so' % (libname, platform.system().lower()), recursive=True)
    if required and not libraries:
        print('error: could not find %s library that matches the current python version' % libname)
        sys.exit(1)
    return [path.basename(filename) for filename in libraries]

if (INTEGRATE_SAMSEG == 'OFF'):
  setup(
      distclass=BinaryDistribution,
      name='pyfs',
      description='A set of python packages to facilitate tools in the FreeSurfer neuroimaging software',
      author='Laboratory for Computational Neuroimaging',
      packages=packages,
      package_data={
          'fsbindings': find_libs('fsbindings') + find_libs('labelfusion', required=False),
          'gems': find_libs('gemsbindings'),
      },
      install_requires=requirements,
      include_package_data=True
   )
else:
  setup(
      distclass=BinaryDistribution,
      name='pyfs',
      description='A set of python packages to facilitate tools in the FreeSurfer neuroimaging software',
      author='Laboratory for Computational Neuroimaging',
      packages=packages,
      package_data={
          'fsbindings': find_libs('fsbindings') + find_libs('labelfusion', required=False),
      },
      install_requires=requirements,
      include_package_data=True
   )

