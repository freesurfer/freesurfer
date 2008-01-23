#!/bin/tcsh -f

set ID='$Id: create_targz.csh,v 1.20 2008/01/23 09:35:01 nicks Exp $'

unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

setenv PLATFORM     $1
setenv RELEASE_TYPE $2
setenv HOSTNAME `hostname -s`

setenv SPACE_FS /space/freesurfer
setenv LOCAL_FS /usr/local/freesurfer
# if /space/freesurfer is down, or if there is a need to run
# the build scripts outside of /usr/local/freesurfer, then 
# the var USE_SPACE_MINERVA can be set
if ($?USE_SPACE_MINERVA) then
  setenv SPACE_FS /space/minerva/1/users/nicks
  setenv LOCAL_FS /space/minerva/1/users/nicks/build/install/${HOSTNAME}
endif

if ("$PLATFORM" == "rh7.3") then
    if ("${HOSTNAME}" != "martinos01" ) then
        echo "must run on machine martinos01"
        exit 1
    endif
else if ("$PLATFORM" == "rh9") then
    if ("${HOSTNAME}" != "kani" ) then
        echo "must run on machine kani"
        exit 1
    endif
else if ("$PLATFORM" == "centos4") then
    if ("${HOSTNAME}" != "fishie" ) then
        echo "must run on machine fishie"
        exit 1
    endif
else if ("$PLATFORM" == "centos4_x86_64") then
    if ("${HOSTNAME}" != "minerva" ) then
        echo "must run on machine minerva"
        exit 1
    endif
else if ("$PLATFORM" == "leopard") then
    if ("${HOSTNAME}" != "storm" ) then
        echo "must run on machine storm"
        exit 1
    endif
else if ("$PLATFORM" == "leopard-i686") then
    if ("${HOSTNAME}" != "sleet" ) then
        echo "must run on machine sleet"
        exit 1
    endif
else if ("$PLATFORM" == "tiger") then
    if ("${HOSTNAME}" != "rain" ) then
        echo "must run on machine rain"
        exit 1
    endif
else
    echo "Usage:"
    echo "$0 <platform> <release_type>"
    echo "where <platform> is rh7.3, rh9, centos4, centos4_x86_64, leopard, leopard-i686 or tiger"
    echo "and <release_type> is either dev, or stable-pub"
    exit 1
endif

if (("$RELEASE_TYPE" != "dev") && ("$RELEASE_TYPE" != "stable-pub") && ("$RELEASE_TYPE" != "stable3-pub")) then
  echo "ERROR: release_type is either dev or stable-pub"
  exit 1
endif

echo "cd ${LOCAL_FS}"
cd ${LOCAL_FS}

if ( ("$PLATFORM" == "leopard") || ("$PLATFORM" == "leopard-i686") || ( "$PLATFORM" == "tiger") ) then
  if (-e /Users/Shared/tmp/$RELEASE_TYPE) \
    rm -Rf /Users/Shared/tmp/$RELEASE_TYPE
  echo "cp -r $RELEASE_TYPE /Users/Shared/tmp"
  cp -r $RELEASE_TYPE /Users/Shared/tmp
  echo "cd /Users/Shared/tmp"
  cd /Users/Shared/tmp
endif

if (-e freesurfer) rm freesurfer
set cmd=(ln -s $RELEASE_TYPE freesurfer)
echo $cmd
$cmd

setenv BUILD_STAMP "`cat ${LOCAL_FS}/${RELEASE_TYPE}/build-stamp.txt`"
setenv FILENAME ${BUILD_STAMP}-full
setenv TARNAME ${FILENAME}.tar
echo creating $TARNAME...
# the -h flag passed to tar is critical!  
# it follows the links to the libraries.
tar -X ${SPACE_FS}/build/scripts/exclude_from_targz -hcvf $TARNAME freesurfer
echo gzipping $TARNAME...
gzip $TARNAME
set cmd=(mv $TARNAME.gz ${SPACE_FS}/build/pub-releases)
echo $cmd
$cmd
set cmd=(chmod g+w ${SPACE_FS}/build/pub-releases/$TARNAME.gz)
echo $cmd
$cmd
if (-e /Users/Shared/tmp/$RELEASE_TYPE) rm -Rf /Users/Shared/tmp/$RELEASE_TYPE

# cleanup useless link
cd ${LOCAL_FS}
if (-e freesurfer) rm freesurfer

# for the Mac, create_dmg.csh creates the .dmg file from the .tar.gz file.
# create_dmg must be executed separately, since root access is 
# required to chown the files to root.
# for Linux, the script postprocess_targz.csh does the same (and also
# requires root access).

