-------------------------------------------------------------------
	Talairaching - Automatic Failure Detection
-------------------------------------------------------------------


METHOD
------

One of the first steps in the processing of anatomical data consists 
in the alignment of the input images with the Talairach atlas. 
A 3 by 4 transform matrix is applied to the input brain, the 12 parameters 
indicating the shifts, rotations, scaling and shears of a linear alignment. 
The severe misalignment of the input brain resulting from a failure in the 
Talairaching is manifest by values of the transform parameters T that are 
outside of the normal range of variation.
In order to automatically detect a wrong Talairach transform for a given 
subject S, we use the 9 components of the transform matrix (parameter vector T) 
corresponding to the rotation and scaling (the 3 translation parameters 
are not taken into account) to compute the probability P(T | mu, sigma) 
of this Talairach transform, mu and sigma being respectively the mean and 
covariance matrix of all the parameter vectors Ti computed for the training set.   

Given the distribution of the probabilities of the Talairach transform parameter 
vectors for the training set, a p-value can also be computed for the subject S. 
(For further reference see dev/docs/Automatic_Failure_Detection.doc)  



HOW TO USE “talairaching_afd.m”
-------------------------------

First, the FREESURFER_HOME variable should be set.

Usage:  [proba, p-value] = talairaching_afd( subjectname, threshold, <statsDirectory>);

The input parameters are:  	
- the subject’s file name
- the p-value threshold below which the Talairach transform is likely to be wrong. 
  This threshold can be adjusted by the user depending on the level of correlation needed.
- the directory where the statistics tables can be found. These table are required to 
  compute the probability and p-value. The training set used by default if statsDirectory 
  is not specified is the buckner data set.

The function returns the probability and the p-value computed for the subject’s 
Talairach transform. If the computed p-value is below the threshold set by the user, 
the Talairach transform is detected as wrong.

(Note: the program looks for the subject’s talairach transform matrix “talairach.xfm”, 
which should be located in subjectname/mri/transform/, if this file doesn’t exist the 
probability cannot be computed.)

exemple:  >> [p, pval] = talairaching_afd('SubjectsDir/bert1', 0.01); 
	  Talairach Transform: OK (p=0.56, pval=0.211)


HOW TO USE "Talairaching_dir_afd.m"
----------------------------------


Usage:  [probas] = talairaching_dir_afd( subjectsDir, threshold, <statsDirectory>);

This function works as the previous one, expect that it checks every subject in the 
directory "subjectsDir".

exemple: >> P=talairaching_dir_afd( 'subjectsDir', 0.01);
            Talairach Transform: bert1 OK, (p=0.56, pval=0.211)
	    Talairach Transform: bert2 failed,  (p=0.14, pval=0.00507)
		...


 	        

