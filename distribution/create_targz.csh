#!/bin/tcsh -ef

set ID='$Id: create_targz.csh,v 1.10 2006/03/04 00:00:26 nicks Exp $'

unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

setenv PLATFORM     $1
setenv RELEASE_TYPE $2
setenv SPACE_FREESURFER /space/freesurfer

if ("$PLATFORM" == "rh7.3") then
    if ("`uname -n`" != "martinos01" ) then
        echo "must run on machine martinos01"
        exit 1
    endif
else if ("$PLATFORM" == "rh9") then
    if ("`uname -n`" != "kani" ) then
        echo "must run on machine kani"
        exit 1
    endif
else if ("$PLATFORM" == "centos4") then
    if ("`uname -n`" != "fishie" ) then
        echo "must run on machine fishie"
        exit 1
    endif
else if ("$PLATFORM" == "centos4_x86_64") then
    if ("`uname -n`" != "minerva" ) then
        echo "must run on machine minerva"
        exit 1
    endif
else if ("$PLATFORM" == "tiger") then
    if ("`uname -n`" != "storm.nmr.mgh.harvard.edu" ) then
        echo "must run on machine storm"
        exit 1
    endif
else
    echo "Usage:"
    echo "$0 <platform> <release_type>"
    echo "where <platform> is rh7.3, rh9, centos4, centos4_x86_64, or tiger"
    echo "and <release_type> is either dev, or stable-pub"
    exit 1
endif

if (( "$RELEASE_TYPE" != "dev") && ( "$RELEASE_TYPE" != "stable-pub")) then
  echo "ERROR: release_type is either dev or stable-pub"
  exit 1
endif

if ( "$PLATFORM" == "tiger") then
  if (-e /Users/Shared/tmp/$RELEASE_TYPE) \
    rm -Rf /Users/Shared/tmp/$RELEASE_TYPE
  cp -R $RELEASE_TYPE /Users/Shared/tmp
  cd /Users/Shared/tmp
endif

cd /usr/local/freesurfer
if (-e freesurfer) rm freesurfer
set cmd=(ln -s $RELEASE_TYPE freesurfer)
echo $cmd
$cmd

setenv BUILD_STAMP "`cat /usr/local/freesurfer/${RELEASE_TYPE}/build-stamp.txt`"
setenv FILENAME ${BUILD_STAMP}-full
setenv TARNAME ${FILENAME}.tar
echo creating $TARNAME...
# the -h flag passed to tar is critical!  it follows the links to the libraries.
tar -X ${SPACE_FREESURFER}/build/scripts/exclude_from_targz -hcvf $TARNAME freesurfer
echo gzipping $TARNAME...
gzip $TARNAME
set cmd=(mv $TARNAME.gz ${SPACE_FREESURFER}/build/pub-releases)
echo $cmd
$cmd
set cmd=(chmod g+w ${SPACE_FREESURFER}/build/pub-releases/$TARNAME.gz)
echo $cmd
$cmd
if (-e /Users/Shared/tmp/$RELEASE_TYPE) rm -Rf /Users/Shared/tmp/$RELEASE_TYPE

# cleanup useless link
cd /usr/local/freesurfer
if (-e freesurfer) rm freesurfer

# for the Mac, create_dmg.csh creates the .dmg file from the .tar.gz file.
# create_dmg must be executed separately, since root access is 
# required to chown the files to root.
# for Linux, the script postprocess_targz.csh does the same (and also
# requires root access).

