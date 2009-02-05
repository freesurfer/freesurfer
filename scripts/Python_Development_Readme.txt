$Id: Python_Development_Readme.txt,v 1.1 2009/02/05 21:49:02 krish Exp $

Introduction
This document is to be used for developers using Python to write scripts and/or extend Freesurfer.

Details 
Usually a script consists of a number of modules and functions and part of the script can be written as an utility function. The main benefit is other scripts in future can use these scripts. For example, currently the 3 scripts
 - datastruct_utils.py
 - subject_info.py
 - misc.py 
exist. The functions inside them are extensively used by a*stats2table scripts. 

When you need to add a new function to any of the above python files ( or a new file), write the function or the class and add the name to __all__ class on top of the file. ( This ensures that the function/class will be included when the script which imports them does so) 

When you need to add a new Python file, create the file and add an __all__ statement ( similar to the above 3 scripts ) which inform what functions to export. Suppose you create foo.py, edit fsutils.py and do a 

from foo import *

so that the script which you'll eventually write just imports fsutils

import fsutils

Rationale
This whole thing reduces namespace cluttering. To the script only one module fsutils.py exist, but you have a few files which each serve a different purpose.
