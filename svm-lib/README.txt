Created by: Polina Golland polina@mit.edu
Created on: 4/4/5


This file explains how to use the classification scripts built around
the SVM code. Following the steps below, you will be able to train a
classifier, visualize its gradient and perform jackknife-style
cross-validation.  It should be self-explanatory, and if it's not,
please let me know. I am happy to change and add comments to help
people to use these scripts.

1. Make sure you have environment variable SVM_DIR defined and
pointing to where the software is located in your system. If you
didn't change the directory structure of the distribution package,
it's one directory up from this README.txt file. The easiest way to do
this is to include it in your .cshrc or .tcshrc file. For example, 

setenv SVM_DIR /afs/csail.mit.edu/u/p/polina/Software/SVM-BIRN


2. Create a working directory. Copy
$SVM_DIR/examples/classification_setup.example to that directory and
edit it to describe your study. 

Note: you will probably have to source subjects.csh or some other file
that defines your groups, where the files are stored, etc. The
definitions within classification_setup.example can certainly involve
your other environment variables that specify the details of the
study.


3. Cd into your working directory and type

source classification_setup.example

The script will run and print a bunch of diagnostics. If it runs into
troubles (not finding files, not being able to read the data, etc.),
it will hopefully give you an informative error message that you can
use to fix your setup. If it completes successfully, the output
directory will contain the following files:

jackknife_error.m - a text file that contains the results of the
jackknife cross-validation. You can either load it in matlab or look
at it with the more command. Each row starts with the size of the
training set from each group used in jackknifing. The rest of the data
was used to test. The next number is the total jackknife error,
followed by group-specific jackknife errors. So plotting the second
column vs. the first one will tell you if you are pretty much
saturated in terms of improving the accuracy, or if the error is still
decreasing rapidly as you add more examples.

study_weights - a file of the same format as your raw data files,
containing the gradient of the classifier, i.e., the vector along
which you should change the examples in the second group to make them
look more like the first group.


The script will also create a bunch of files in the files
directory. Those are temporary files that are useful for
debugging. You can pretty much ignore them (remove them if you want),
if you are not trying to reconstruct what went wrong.


