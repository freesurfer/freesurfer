Survival
========

SStat (survival stats) matlab toolbox for Survival (time-to-event) analysis. Implements both univariate and mass-univariate analyses. Jorge Luis Bernal Rusiel, 2013.

If you use these tools in your analyses please cite:

Sabuncu M.R., Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., 2014. Event time analysis of longitudinal neuroimage data, Neuroimage doi: 10.1016/j.neuroimage.2014.04.015.

These Matlab tools are freely distributed and intended to help neuroimaging researchers when analyzing time-to-event neuroimage data with time-independent and perhaps time-dependent (longitudinal) covariates as predictors of the event timing (eg. disease onset). They provide functionality for Kaplan-Meier curves estimation and plot, proportional hazards(PH) test, hazard ratio(HR) calculation, model specification, parameter estimation and inference using Cox regression models. They are specially targeted to be used with Freesurfer's data but can be used with any other data as long as they are loaded into Matlab and put in the appropriate format. 

Take a look at the lme toolbox's wiki to see how the Neuroimage data (eg. cortical thickness measurements along the cortical mantle) can be generated with Freesurfer. Also read thoroughly the paper referenced above to understand the theory behind the joint lme-Cox statistical analysis. A practical example using ADNI data for advanced statistical analyses of longitudinal neuroimage predictors of conversion from mild cognitive impairment to Alzheimer's disease is also presented.

The main reference (highly recommended) used for the implementation of the Cox models was:

Kleinbaum, D.G. and Klein, M., 2012. Survival Analysis. Springer.

Please read the header of each script to understand its functionality as well as its inputs and outputs. Matlab files with some example data sets ready for the analysis with this toolbox can be found within the ex_data directory. For example you can load the file vet.mat into Matlab and perform a univariate analysis with the Cox PH regression model:

[stats,st] = SStat_CoxPH(X,d,t); 

When combining lme and Cox models (joint analysis) the subject specific random effects estimated from the longitudinal repeated measures can be computed after a typical lme analysis with:

lme_mass_rfx

then the extended Cox design matrix can be obtained using:

SStat_lmeX_ext

You can also post questions about this toolbox to Freesurfer mailing list at freesurfer[at]nmr.mgh.harvard.edu.





