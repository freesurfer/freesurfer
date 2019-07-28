import os
import sys
from setuptools import setup
from glob import glob


setup(
    name='sscnn_skullstripping',
    packages=['sscnn_skullstripping', 'sscnn_skullstripping.predict','sscnn_skullstripping.deeplearn_utils' ],
    url='https://surfer.nmr.mgh.harvard.edu/',
    license='FreeSurfer Software License Agreement',
    author='Amod S. Jog',
    author_email='amoddoma@gmail.com',
    description='SSCNN: Fast Multi-view 2D CNN Infant MRI Skullstripping',
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Operating System :: MacOS :: MacOS X',
                 'Operating System :: POSIX :: Linux',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Topic :: Scientific/Engineering'],
    install_requires=['keras',
                      'tensorflow',
                      'numpy',
                      'scipy',
                      'nibabel',
                      'matplotlib',
                      'seaborn',
                      'configparser'],
    dependency_links=[],
    provides=['sscnn_skullstripping'],
    package_data={'sscnn_skullstripping': ['model_files/*.h5']},
    scripts=glob('bin/*')

)
