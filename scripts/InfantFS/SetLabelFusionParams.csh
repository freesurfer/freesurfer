#!/bin/tcsh -ef    

#setenv CMDDIR `pwd`
setenv CMDDIR $FREESURFER_HOME/scripts/LabelFusionBetterDist/
setenv MLAB /autofs/cluster/matlab/8.0/
setenv MRFparams "1 5"

setenv LFOPTIONS1 "2 41 3 42" # new -- cerebral WM hack
setenv LFOPTIONS2 "7 46 8 47" # new -- cerebellar WM hack
setenv LFOPTIONS3 "12 51 13 52" # new -- basal ganglia hack

setenv NEGLFOPTIONS1 "-2 -41 -3 -42" # new
setenv NEGLFOPTIONS2 "-7 -46 -8 -47" # new
setenv NEGLFOPTIONS3 "-12 -51 -13 -52" # new

# set LFOPTIONS = (2 3 41 42) # new
# set LFOPTIONS = (41 2 42 3)

setenv MAXLAB 3
setenv BIASFIELDORDER 4
setenv BETA 0.3 # 0.1
