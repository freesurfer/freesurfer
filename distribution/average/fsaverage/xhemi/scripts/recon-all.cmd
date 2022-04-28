

#---------------------------------
# New invocation of recon-all Fri Nov 13 11:35:23 EST 2015 
#--------------------------------------------
#@# Talairach Fri Nov 13 11:35:32 EST 2015

 mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz --n 1 --proto-iters 1000 --distance 200 --stop 1e-4 --shrink 2 


 talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm 

talairach_avi log file is transforms/talairach_avi.log...

 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

