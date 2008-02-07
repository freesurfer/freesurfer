Created by: Polina Golland polina@mit.edu
Created on: 4/4/5

This file explains how to install the classification scripts built
around the SVM code. Following the steps below, you will be able to
install the software and then use it as described in README.txt It
should be self-explanatory, and if it's not, please let me know. I am
happy to change and add comments to help people to use these scripts.

1. Copy the entire SvmClass directory to where you typically place
code. 

2. Compile the library and the applications.

type 

make clean 
make depend
make

in each of these two directories:

SvmClass/code/dist/svm-lib/src
SvmClass/code/dist/app/src

3. Once you complied the code, it is ready to be used.

SvmClass/examples

contains an example script.

