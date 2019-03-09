import os
import sys
from setuptools import setup
from glob import glob


setup(
    name='psacnn_brain_segmentation',
    packages=['psacnn_brain_segmentation', 'psacnn_brain_segmentation.predict','psacnn_brain_segmentation.deeplearn_utils' ],
    url='https://surfer.nmr.mgh.harvard.edu/',
    license='FreeSurfer Software License Agreement',
    author='Amod S. Jog',
    author_email='amoddoma@gmail.com',
    description='PSACNN: Pulse Sequence Adaptive CNN for Whole Brain Segmentation',
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
    provides=['psacnn_brain_segmentation'],
    package_data={'psacnn_brain_segmentation': ['model_files/*.h5']},
    scripts=glob('bin/*')

)
