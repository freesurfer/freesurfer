#-----------------------------------------
#@# Cortical Parc 2 lh Thu Dec 30 17:40:05 EST 2010

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 fsaverage lh ../surf/lh.sphere.reg /usr/local/freesurfer/stable5/average/lh.destrieux.simple.2009-07-29.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Cortical Parc 2 rh Thu Dec 30 17:40:58 EST 2010

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.mgz -seed 1234 fsaverage rh ../surf/rh.sphere.reg /usr/local/freesurfer/stable5/average/rh.destrieux.simple.2009-07-29.gcs ../label/rh.aparc.a2009s.annot 

