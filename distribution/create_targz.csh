#!/bin/tcsh -f

set ID='$Id: create_targz.csh,v 1.34 2013/05/22 22:04:42 nicks Exp $'

unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

setenv PLATFORM     $1
setenv RELEASE_TYPE $2
setenv HOSTNAME `hostname -s`

setenv OSTYPE `uname -s`
if ("$OSTYPE" == "linux") setenv OSTYPE Linux
if ("$OSTYPE" == "Linux") setenv OSTYPE Linux
if ("$OSTYPE" == "darwin") setenv OSTYPE Darwin
if ("$OSTYPE" == "Darwin") setenv OSTYPE Darwin

setenv SPACE_FS /space/freesurfer
setenv LOCAL_FS /usr/local/freesurfer
# if /space/freesurfer is down, or if there is a need to install
# outside of /usr/local/freesurfer, 
# then these can override
if ($?USE_SPACE_FS) then
  setenv SPACE_FS $USE_SPACE_FS
endif
if ($?USE_LOCAL_FS) then
  setenv LOCAL_FS $USE_LOCAL_FS
endif
# if /space/freesurfer is down, or if there is a need to run
# the build scripts outside of /usr/local/freesurfer, then 
# the var USE_SPACE_MINERVA can be set
if ($?USE_SPACE_MINERVA) then
  setenv SPACE_FS /space/minerva/1/users/nicks
  setenv LOCAL_FS /space/minerva/1/users/nicks/build/install/${HOSTNAME}
endif

if ("$PLATFORM" == "centos6_x86_64") then
    if ("${HOSTNAME}" != "terrier" ) then
        echo "must run on machine terrier"
        exit 1
    endif
else if ("$PLATFORM" == "centos5_x86_64") then
    if ("${HOSTNAME}" != "swan" ) then
        echo "must run on machine swan"
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
else if ("$PLATFORM" == "leopard-ppc") then
    if ("${HOSTNAME}" != "storm" ) then
        echo "must run on machine storm"
        exit 1
    endif
else if ("$PLATFORM" == "snow_leopard") then
    if ("${HOSTNAME}" != "sleet" ) then
        echo "must run on machine sleet"
        exit 1
    endif
else if ("$PLATFORM" == "snowleopard-i686") then
    if ("${HOSTNAME}" != "sleet" ) then
        echo "must run on machine sleet"
        exit 1
    endif
else if ("$PLATFORM" == "leopard-i686") then
    if ("${HOSTNAME}" != "hima" ) then
        echo "must run on machine hima"
        exit 1
    endif
else if ("$PLATFORM" == "mountain_lion") then
    if ("${HOSTNAME}" != "hima" ) then
        echo "must run on machine hima"
        exit 1
    endif
else if ("$PLATFORM" == "lion") then
    if ("${HOSTNAME}" != "gust" ) then
        echo "must run on machine gust"
        exit 1
    endif
else if ("$PLATFORM" == "tiger-ppc") then
    if ("${HOSTNAME}" != "mist" ) then
        echo "must run on machine mist"
        exit 1
    endif
else if ("$PLATFORM" == "tiger-i686") then
    if ("${HOSTNAME}" != "sleet" ) then
        echo "must run on machine sleet"
        exit 1
    endif
else
    echo "Usage:"
    echo "$0 <platform> <release_type>"
    echo "where <platform> is centos6_x86_64, centos4, centos4_x86_64,"
    echo "snow_leopard, lion, or mountain_lion"
    echo "and <release_type> is either dev, or stable-pub"
    exit 1
endif

if (("$RELEASE_TYPE" != "dev") && \
    ("$RELEASE_TYPE" != "stable-pub") \
    && ("$RELEASE_TYPE" != "stable3-pub")) then
  echo "ERROR: release_type is either dev or stable-pub"
  exit 1
endif

echo "cd ${LOCAL_FS}"
cd ${LOCAL_FS}

if ( ("$PLATFORM" == "leopard-ppc") || \
    ("$PLATFORM" == "leopard-i686") || \
    ( "$PLATFORM" == "tiger-ppc") ) then
  if (-e /Users/Shared/tmp/$RELEASE_TYPE) \
    rm -Rf /Users/Shared/tmp/$RELEASE_TYPE
  echo "cp -r $RELEASE_TYPE /Users/Shared/tmp"
  cp -r $RELEASE_TYPE /Users/Shared/tmp
  echo "cd /Users/Shared/tmp"
  cd /Users/Shared/tmp
endif

rm -f freesurfer
set cmd=(ln -s $RELEASE_TYPE freesurfer)
echo $cmd
$cmd

if ("$RELEASE_TYPE" == "dev") then
  setenv FILENAME freesurfer-${OSTYPE}-${PLATFORM}-dev
else
  setenv BUILD_STAMP "`cat ${LOCAL_FS}/${RELEASE_TYPE}/build-stamp.txt`"
  setenv FILENAME ${BUILD_STAMP}
endif
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

# check contents of tar file
#${SPACE_FS}/build/scripts/check_manifest.sh ${SPACE_FS}/build/pub-releases/$TARNAME.gz ${SPACE_FS}/build/scripts/manifest/${TARNAME}.gz.manifest

# for the Mac, create_dmg.csh creates the .dmg file from the .tar.gz file.
# create_dmg must be executed separately, since root access is 
# required to chown the files to root.
# for Linux, the script postprocess_targz.csh does the same (and also
# requires root access).

