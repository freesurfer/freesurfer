#!/bin/csh -f

# rebuild_gca_atlas.csh
# This script builds the subcortical atlas from a set of manually labelled 
# training data sets. 
# The user should call the script under a cluster machine (seychelles).
# In addtion, a "subject.csh" script should be available under the 
# $SUBJECTS_DIR/scripts directory, which defines the training subject names
# and the first subject ($ONE_SUBJECT) to be used for building the initial
# template.
# A talairach registration should be generated for this first subject as well,
# in order to align the final atlas to the Talairach space
# The registration file should be put under 
# $SUBJECTS_DIR/$ONE_SUBJECT/mri/transforms/$talairach_manual 
# The script assumes the existence of following files under 
# each subject's mri directory: nu.mgz, brain.mgz, and nu_noneck.mgz
# The final atlas is stored as $SUBJECTS_DIR/average/RB_all_`date +%F`.gca
# "nu_noneck.mgz" is needed for every subject to build the gca with skull

set VERSION = '$Id: rebuild_gca_atlas.csh,v 1.1 2006/09/27 17:49:33 nicks Exp $';

source $SUBJECTS_DIR/scripts/subjects_old.csh

#programs needed
setenv emreg mri_em_register
setenv careg mri_ca_register
setenv canorm mri_ca_normalize
setenv train mri_ca_train

#parameters may need be modified 
set SEG_VOL = seg_edited6.mgz #filename for manual segmentation
set ORIG_VOL = nu.mgz
set maskvol =  brain.mgz #filename for brain mask
set T1_NONECK = nu_noneck.mgz #file to build the atlas gca_with_skull
set talairach_manual = RB_AVERAGE3new.xfm 
#talairach_manual.lta

mkdir -p $SUBJECTS_DIR/average

setenv GCA $SUBJECTS_DIR/average/RB_all_`date +%F`.gca
setenv GCA_ONE $SUBJECTS_DIR/average/RB_one_`date +%F`.gca
setenv GCAwithskull $SUBJECTS_DIR/average/RB_all_withskull_`date +%F`.gca

setenv LTA_ONE talairach_one.lta 
setenv M3D_ONE talairach_one.m3z

setenv LTA talairach.lta
setenv M3D talairach.m3z

echo START GCA TRAINING

set T1_VOL = norm.mgz

#normalize brains
echo INITIAL NORMALIZE
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
        if(-e $mridir/${T1_VOL}) then
           rm $mridir/${T1_VOL}
        endif
        pbsubmit -l nodes=1:new -c "$canorm -mask $mridir/$maskvol -seg $mridir/$SEG_VOL $mridir/$ORIG_VOL noatlas noxform $mridir/${T1_VOL}"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    setenv INNORM $mridir/${T1_VOL}
                set TEST = 1
                echo $subject
                while($TEST)
                                if(-e $INNORM) then
                                                set TEST = 0
                                endif
                                sleep 30
                end
end


##TRAIN FROM ONE_SUBJECT
echo TRAIN_ONE

if( -e $GCA_ONE ) then
  rm $GCA_ONE
endif

if(-e $SUBJECTS_DIR/$ONE_SUBJECT/mri/transforms/${talairach_manual}) then 
pbsubmit -l nodes=1:new -c "$train -prior_spacing 2 -node_spacing 8 -mask $maskvol -parc_dir $SEG_VOL -xform ${talairach_manual}  -T1 ${T1_VOL} $ONE_SUBJECT $GCA_ONE"
else
echo "Initial Talairach Registration Unavailable"
pbsubmit -l nodes=1:new -c "$train -prior_spacing 2 -node_spacing 8 -mask $maskvol -parc_dir $SEG_VOL -T1 ${T1_VOL} $ONE_SUBJECT $GCA_ONE"
endif

#waiting 
echo WAITING
set TEST = 1 
while($TEST) 
	 if(-e $GCA_ONE) then
	      set TEST = 0
	 endif
	 sleep 30
end

# EM_registration by GCA_ONE
echo EM_registration by GCA_ONE
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    mkdir -p $mridir/transforms
    setenv tdir $mridir/transforms
    if ( -e $tdir/$LTA_ONE ) then
	rm $tdir/$LTA_ONE
    endif	
    pbsubmit -l nodes=1:new -c "$emreg -mask $mridir/$maskvol $mridir/$ORIG_VOL $GCA_ONE $mridir/transforms/$LTA_ONE"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv INREG $SUBJECTS_DIR/$subject/mri/transforms/$LTA_ONE
	  set TEST = 1 
	  echo $subject
	  while($TEST) 
				if(-e $INREG) then
					   set TEST = 0
				endif
				sleep 30
		end
end

# normalization by GCA_ONE
echo Normalization by GCA_ONE
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    if ( -e $mridir/$T1_VOL ) then
	rm $mridir/$T1_VOL
    endif
  pbsubmit -l nodes=1:new -c "$canorm -mask $mridir/$maskvol $mridir/$ORIG_VOL $GCA_ONE $mridir/transforms/$LTA_ONE $mridir/$T1_VOL"
end

#wait untill the normalization is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv INNORM $SUBJECTS_DIR/$subject/mri/$T1_VOL
		set TEST = 1 
		echo $subject
		while($TEST) 
				if(-e $INNORM) then
						set TEST = 0
				endif
				sleep 30
		end
end



# CA_registration by GCA_ONE
echo CA_registration by GCA_ONE
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    if ( -e  $mridir/transforms/$M3D_ONE ) then
	rm $mridir/transforms/$M3D_ONE
    endif 	
    pbsubmit -l nodes=1:new -c "$careg -smooth 1.0 -levels 2 -mask $mridir/$maskvol -T $mridir/transforms/$LTA_ONE $mridir/$T1_VOL $GCA_ONE $mridir/transforms/$M3D_ONE"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv INREG $SUBJECTS_DIR/$subject/mri/transforms/$M3D_ONE
	  set TEST = 1 
	  echo $subject
	  while($TEST) 
				if(-e $INREG) then
					   set TEST = 0
				endif
				sleep 30
		end
end

##TRAIN FROM SEGMENTED_SUBJECTS USING M3D_ONE
echo TRAIN_ALL
if (-e $GCA ) then
    rm -f $GCA
endif 

pbsubmit -l nodes=1:new -c "$train -prior_spacing 2 -node_spacing 4 -mask $maskvol -parc_dir $SEG_VOL -xform $M3D_ONE -T1 $T1_VOL $SUBJECTS $GCA"

#waiting 
echo WAITING
set TEST = 1 
while($TEST) 
	 if(-e $GCA) then
	    set TEST = 0
	 endif
	 sleep 30
end

#REGISTER ALL BRAIN TO GCA
# EM_registration by GCA
echo FIRST EM registration by GCA
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    setenv tdir $mridir/transforms
    if ( -e $tdir/$LTA ) then
	rm $tdir/$LTA
    endif 	
    pbsubmit -l nodes=1:new -c "$emreg -mask $mridir/$maskvol $mridir/$ORIG_VOL $GCA $mridir/transforms/$LTA"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv INREG $SUBJECTS_DIR/$subject/mri/transforms/$LTA
	  set TEST = 1 
	  echo $subject
	  while($TEST) 
				if(-e $INREG) then
					   set TEST = 0
				endif
				sleep 30
		end
end

# normalization by GCA
echo FIRST Normalization by GCA
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    if ( -e $mridir/$T1_VOL ) then
	rm $mridir/$T1_VOL
    endif 	
  pbsubmit -l nodes=1:new -c "$canorm -mask $mridir/$maskvol $mridir/$ORIG_VOL $GCA $mridir/transforms/$LTA $mridir/$T1_VOL"
end

#wait untill the normalization is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv INNORM $SUBJECTS_DIR/$subject/mri/$T1_VOL
		set TEST = 1 
		echo $subject
		while($TEST) 
				if(-e $INNORM) then
						set TEST = 0
				endif
				sleep 30
		end
end



# CA_registration by GCA
echo FIRT CA_registration by GCA
foreach subject ($SUBJECTS)
    setenv mridir $SUBJECTS_DIR/$subject/mri
    if ( -e $mridir/transforms/$M3D ) then
	rm $mridir/transforms/$M3D
    endif
	
    pbsubmit -l nodes=1:new -c "$careg -smooth 1.0 -mask $mridir/$maskvol -T $mridir/transforms/$LTA $mridir/$T1_VOL  $GCA $mridir/transforms/$M3D"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
    setenv INREG $SUBJECTS_DIR/$subject/mri/transforms/$M3D
	  set TEST = 1 
	  echo $subject
	  while($TEST) 
				if(-e $INREG) then
					   set TEST = 0
				endif
				sleep 30
		end
end

#RETRAIN GCA USING M3D
echo TRAIN_ALL_SECOND
rm -f $GCA

pbsubmit -l nodes=1:new -c "$train -prior_spacing 2 -node_spacing 4 -mask  $maskvol -parc_dir $SEG_VOL -xform $M3D -T1 $T1_VOL $SUBJECTS $GCA"

#waiting 
echo WAITING
set TEST = 1 
while($TEST) 
	 if(-e $GCA) then
	    set TEST = 0
	 endif
	 sleep 30
end

rm -f $GCAwithskull

pbsubmit -l nodes=1:new -c "$train -prior_spacing 2 -node_spacing 4 -parc_dir $SEG_VOL -xform $LTA -T1 $T1_NONECK $SUBJECTS $GCAwithskull"

#waiting 
echo WAITING
set TEST = 1 
while($TEST) 
	 if(-e $GCA) then
	    set TEST = 0
	 endif
	 sleep 30
end

echo DONE
