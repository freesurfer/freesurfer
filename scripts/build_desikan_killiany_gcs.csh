#!/bin/tcsh -ef

# For building a new morph-based cortical parcellation 
# atlas for each hemisphere

# The official SUBJECTS_DIR for the training subjects is:
# /autofs/space/amaebi_026/users/buckner_cortical_atlas
# SUBJECTS_DIR should be set to that (or some other training set) 
# prior to running this script.  The .gcs file is written to
# the 'average' directory of SUBJECTS_DIR, along with a log file.
# The script 'scripts/subjects.csh' in SUBJECTS_DIR must set
# the var SUBJECTS to the list of training subjects.

set HEMI = $1
if ("x$HEMI" != "xrh" && "x$HEMI" != "xlh") then
    echo "Usage: $0 <hemi>"
    echo "where <hemi> is rh or lh"
    exit 1
endif
set DATE = "`date '+%Y-%m-%d'`"
set GCS = ${HEMI}.curvature.buckner40.filled.desikan_killiany.${DATE}.gcs
set LOG = $SUBJECTS_DIR/average/mris_ca_train-${HEMI}-${DATE}.log
source scripts/subjects.csh
set cmd=(mris_ca_train \
  -t $FREESURFER_HOME/average/colortable_desikan_killiany.txt \
  ${HEMI} \
  sphere.reg \
  aparc_edited \
  $SUBJECTS \
  $SUBJECTS_DIR/average/${GCS})

mris_ca_train --all-info >& $LOG
echo "" >> $LOG
echo "Start: `date`" >> $LOG
echo "" >> $LOG
echo "${cmd} \n" >> $LOG

$cmd >>& $LOG

echo "" >> $LOG
echo "End:   `date`" >> $LOG
echo "" >> $LOG
