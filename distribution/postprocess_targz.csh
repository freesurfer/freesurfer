#!/bin/tcsh -f

set echo=1
cd $MINERVA_HOME/tmp
if (-e freesurfer) sudo rm -Rf freesurfer
sudo tar zxvf /space/freesurfer/build/pub-releases/$1.tar.gz
if ($status) exit 1
sudo chmod -R a+rw freesurfer
if ($status) exit 1
sudo chmod -R a-w freesurfer/subjects/fsaverage
if ($status) exit 1
sudo chown -R root freesurfer
if ($status) exit 1
sudo chown -R root freesurfer/*.*
if ($status) exit 1
sudo chgrp -R root freesurfer
if ($status) exit 1
sudo chgrp -R root freesurfer/*.*
if ($status) exit 1
sudo strip freesurfer/bin/* >& /dev/null
# don't check strip return status, as strip returns error because 
# it finds files that aren't binaries, but thats ok
tar cvf $1.tar freesurfer
if ($status) exit 1
gzip $1.tar
if ($status) exit 1
md5sum $1.tar.gz >> /space/freesurfer/build/pub-releases/md5sum.txt
mv $1.tar.gz /space/freesurfer/build/pub-releases/
if ($status) exit 1
exit 0
