#!/bin/csh -f

#For building a new surface morph for one hemisphere

set hemi = $1

setenv SUBJECTS_DIR /autofs/space/birn_044/users/buckner_surf_atlas_rebuild
source $SUBJECTS_DIR/scripts/subjects.csh

setenv AVERAGE $SUBJECTS_DIR/average

setenv ONE_SUBJECT $RB_MID[1]
set ONE_SUBJECT = $RB_MID[1]

############################################
set MORPH_SUFFIX = curvature.filled.buckner40
############################################

# make first template
echo FIRST TEMPLATE
if(-e $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif) then
  rm $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif
endif
#is the -norot right in the following? or does it matter? It doesn't!
pbsubmit  -c "mris_make_template  -norot ${hemi} sphere $ONE_SUBJECT $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif"

#waiting 
echo WAITING
set TEST = 1 
setenv INSURF $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif
while($TEST) 
	 if(-e $INSURF) then
			set TEST = 0
	 endif
	 sleep 30
end

#register brains
echo FIRST REGISTRATION
foreach subject ($SUBJECTS)
    	setenv SURF $SUBJECTS_DIR/$subject/surf
	if(-e $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg) then
	   rm $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg
	endif   
	pbsubmit  -c "mris_register -w 0 -curv $SURF/${hemi}.sphere $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
		setenv INSURF $SUBJECTS_DIR/$subject/surf/${hemi}.sphere.${MORPH_SUFFIX}.reg
		set TEST = 1 
		echo $subject
		while($TEST) 
				if(-e $INSURF) then
						set TEST = 0
				endif
				sleep 30
		end
end


#regenerate the template
echo SECOND TEMPLATE
mv $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif $AVERAGE/${hemi}.average.${MORPH_SUFFIX}_${ONE_SUBJECT}.tif
pbsubmit  -c "mris_make_template -norot ${hemi} sphere.${MORPH_SUFFIX}.reg $SUBJECTS $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif"

#waiting 
echo WAITING
set TEST = 1 
setenv INSURF $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif
while($TEST) 
	 if(-e $INSURF) then
			set TEST = 0
	 endif
	 sleep 30
end

# register brains
echo SECOND REGISTRATION
foreach subject ($SUBJECTS)
    setenv SURF $SUBJECTS_DIR/$subject/surf
    rm $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg
    pbsubmit  -c "mris_register -w 0 -curv $SURF/${hemi}.sphere $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
		setenv INSURF $SUBJECTS_DIR/$subject/surf/${hemi}.sphere.${MORPH_SUFFIX}.reg
		set TEST = 1 
		echo $subject
		while($TEST) 
				if(-e $INSURF) then
						set TEST = 0
				endif
				sleep 30
		end
end


#final template
echo FINAL TEMPLATE
mv $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif $AVERAGE/${hemi}.average.${MORPH_SUFFIX}_second.tif
pbsubmit  -c "mris_make_template  -norot  ${hemi} sphere.${MORPH_SUFFIX}.reg $SUBJECTS $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif" 

#waiting 
echo WAITING
set TEST = 1 
setenv INSURF $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif
while($TEST) 
	 if(-e $INSURF) then
			set TEST = 0
	 endif
	 sleep 30
end

#final registration
echo FINAL REGISTRATION
foreach subject ($SUBJECTS)
     setenv SURF $SUBJECTS_DIR/$subject/surf
     rm $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg
     pbsubmit  -c "mris_register -w 0 -curv $SURF/${hemi}.sphere $AVERAGE/${hemi}.average.${MORPH_SUFFIX}.tif $SURF/${hemi}.sphere.${MORPH_SUFFIX}.reg"
end

#wait untill the registration is done
echo WAITING
foreach subject ($SUBJECTS)
		setenv INSURF $SUBJECTS_DIR/$subject/surf/${hemi}.sphere.${MORPH_SUFFIX}.reg
		set TEST = 1 
		echo $subject
		while($TEST) 
				if(-e $INSURF) then
						set TEST = 0
				endif
				sleep 30
		end
end

echo DONE
