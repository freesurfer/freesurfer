#--------------------------------------------
#@# Cortical ribbon mask Thu Jan 12 13:10:41 EST 2012

 mris_volmask --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon ni.sym.b40.i30 

#-----------------------------------------
#@# AParc-to-ASeg Thu Jan 12 13:15:04 EST 2012

 mri_aparc2aseg --s ni.sym.b40.i30 --volmask 


 mri_aparc2aseg --s ni.sym.b40.i30 --volmask --a2009s 

