#!/bin/tcsh -f

set ID='$Id: postprocess_targz.csh,v 1.16 2013/05/22 22:04:42 nicks Exp $'

set echo=1

setenv SPACE_FS /space/freesurfer
# if /space/freesurfer is down, or if there is a need to run
# the build scripts outside of /usr/local/freesurfer, then 
# the var USE_SPACE_MINERVA or USE_SPACE_FS can be set
if ($?USE_SPACE_MINERVA) then
  setenv SPACE_FS /space/minerva/1/users/nicks
endif
if ($?USE_SPACE_FS) then
  setenv SPACE_FS $USE_SPACE_FS
endif

cd /tmp
if ( ! -e scratch ) mkdir scratch
cd scratch

if (-e freesurfer) sudo rm -Rf freesurfer
sudo tar zxvf ${SPACE_FS}/build/pub-releases/$1.tar.gz
if ($status) exit 1

# compress binaries using 'upx', reducing their size even more than strip (3x)
if ($?RUN_UPX) then
  if ( -e /usr/pubsw/bin/upx ) then
    rm -f ELF-files
    foreach f (`ls freesurfer/bin/*`)
      file $f | grep ELF >& /dev/null
      if ( ! $status ) echo $f >> ELF-files
    end
    # UPX seems to munge the _cuda file, causing this error:
    # CUDA Error in file 'devicemanagement.cu'
    # so dont upx those files
    grep -v _cuda ELF-files >> ELF-files-wocuda
    mv ELF-files-wocuda ELF-files
    # now run upx...
    echo "ELF-files:"
    cat ELF-files
    foreach f (`cat ELF-files`)
      sudo /usr/pubsw/bin/upx $f
    end
  endif
endif

sudo chmod -R u+rwX freesurfer
sudo chmod -R go-w freesurfer
sudo chmod -R go+rX freesurfer
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
cd freesurfer/lib/vtk/lib
sudo ${SPACE_FS}/build/scripts/fixvtklibs.csh
cd -
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
