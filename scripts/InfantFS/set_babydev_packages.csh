#!/bin/tcsh -ef

##### scripts location

set FS_HOME = $FREESURFER_HOME
setenv FSSCRIPTSDIR      $FS_HOME/bin/
setenv BABYDEVSCRIPTSDIR $FSSCRIPTSDIR
set BabyDevPath = $BABYDEVSCRIPTSDIR

##### IN ORDER TO ENABLE PICASSO TO RUN

  set FSL_DIR      = /usr/pubsw/packages/fsl/current
  set FSL_DIR2     = /usr/pubsw/packages/fsl/current/bin
  set AFNI_DIR     = /usr/pubsw/packages/AFNI/current
  set N4_DIR       = /usr/pubsw/packages/N4/Build
  set DRAMMS_DIR   = /usr/pubsw/packages/DRAMMS/1.4.4/bin
  set DRAMMS_DIR2  = /usr/pubsw/packages/DRAMMS/1.4.4/lib
  set c3d_DIR      = /usr/pubsw/packages/c3d/bin
  set NiftySeg_DIR = /usr/pubsw/packages/NiftiSeg/build/seg-apps
  set ROBEX_DIR    = /usr/pubsw/packages/ROBEX
  set BSE16a_DIR   = /usr/pubsw/packages/BrainSuite/16a/bin
  set normalizeFOV_DIR = /usr/pubsw/packages/normalizeFOV/1.2/bin
  set picasso_DIR  = /usr/pubsw/packages/picasso/bin

set AllPicassoPath = ()
  foreach dir ( ${FSL_DIR} ${FSL_DIR2} ${AFNI_DIR} ${N4_DIR} ${DRAMMS_DIR} ${DRAMMS_DIR2} ${c3d_DIR} ${NiftySeg_DIR} ${ROBEX_DIR} ${BSE16a_DIR}  ${normalizeFOV_DIR} ${picasso_DIR} )
    if ( -d ${dir} ) then
        set AllPicassoPath = ( ${AllPicassoPath}:${dir} )
        #echo "loaded: ${dir}"
    endif
  end

#####

setenv PATH $PATH${AllPicassoPath}:${BabyDevPath}:${FSSCRIPTSDIR}:$FS_HOME/scripts/infantFS
echo $PATH


#####
##### FOR LABELFUSION
source $BABYDEVSCRIPTSDIR/SetLabelFusionParams.csh

