#!/bin/tcsh -f

set ID='$Id: postprocess_targz.csh,v 1.14 2011/05/03 20:35:09 nicks Exp $'

set echo=1

setenv SPACE_FS /space/freesurfer
# if /space/freesurfer is down, or if there is a need to run
# the build scripts outside of /usr/local/freesurfer, then 
# the var USE_SPACE_MINERVA can be set
if ($?USE_SPACE_MINERVA) then
  setenv SPACE_FS /space/minerva/1/users/nicks
endif

cd /tmp
if ( ! -e scratch ) mkdir scratch
cd scratch

if (-e freesurfer) sudo rm -Rf freesurfer
sudo tar zxvf ${SPACE_FS}/build/pub-releases/$1.tar.gz
if ($status) exit 1
sudo chmod -R u+rw freesurfer
sudo chmod -R go-w freesurfer
if ($status) exit 1
#sudo chmod -R a-w freesurfer/subjects/fsaverage
#if ($status) exit 1
sudo chown -R root freesurfer
if ($status) exit 1
sudo chown -R root freesurfer/*.*
if ($status) exit 1
set ROOTGRP=root
if ("`uname -s`" == "Darwin") set ROOTGRP=wheel
sudo chgrp -R $ROOTGRP freesurfer
if ($status) exit 1
sudo chgrp -R $ROOTGRP freesurfer/*.*
if ($status) exit 1
pwd
tar -X ${SPACE_FS}/build/scripts/exclude_from_targz -cvf $1.tar freesurfer
if ($status) exit 1
gzip $1.tar
if ($status) exit 1
md5sum $1.tar.gz >> ${SPACE_FS}/build/pub-releases/md5sum.txt
sha1sum $1.tar.gz >> ${SPACE_FS}/build/pub-releases/sha1sum.txt
mv $1.tar.gz ${SPACE_FS}/build/pub-releases/
if ($status) exit 1
sudo rm -Rf freesurfer

exit 0
