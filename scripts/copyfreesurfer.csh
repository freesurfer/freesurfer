#! /bin/tcsh -f

#
# copyfreesurfer.csh
#
# REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
#
# Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR
# CVS Revision Info:
#    $Author: nicks $
#    $Date: 2011/03/02 20:16:38 $
#    $Revision: 1.6 $
#
# Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
#
# Terms and conditions for use, reproduction, distribution and contribution
# are found in the 'FreeSurfer Software License Agreement' contained
# in the file 'LICENSE' found in the FreeSurfer distribution, and here:
#
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
#
# Reporting: freesurfer@nmr.mgh.harvard.edu
#
#


set SOURCE_DIR = $argv[1]
set DEST_DIR = $argv[2]

if(! -e $SOURCE_DIR) then
  echo "ERROR: cannot find $SOURCE_DIR"
  exit 1;
endif

set EXCLUDE_FILE = /tmp/exclude-cp
rm -f $EXCLUDE_FILE

pushd $SOURCE_DIR > /dev/null
cd ..

set SOURCE_DIR_REL = `basename $SOURCE_DIR`

find $SOURCE_DIR_REL -name \*.03\* >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*.02\* >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*.01\* >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*~ >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \#\*\# >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*.bu >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*.bak >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*.log >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name \*.old >> $EXCLUDE_FILE
find $SOURCE_DIR_REL -name .xdebug\* >> $EXCLUDE_FILE

rsync --recursive -v --copy-links --exclude-from=$EXCLUDE_FILE $SOURCE_DIR $DEST_DIR


rm -f $EXCLUDE_FILE
