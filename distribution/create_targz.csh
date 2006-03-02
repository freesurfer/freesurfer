#!/bin/tcsh -ef

unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

setenv PLATFORM     $1
setenv RELEASE_TYPE $2

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
else if ("$PLATFORM" == "rhel4") then
    if ("`uname -n`" != "triassic" ) then
        echo "must run on machine triassic"
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
    echo "where <platform> is rh7.3, rh9, rhel4, centos4, centos4_x86_64, or tiger"
    echo "and <release_type> is either dev, or stable-pub"
    exit 1
endif

setenv SPACE_FREESURFER /space/freesurfer
cd $SPACE_FREESURFER
if ( "$PLATFORM" == "centos4") then
    cd centos4.0
else if ( "$PLATFORM" == "centos4_x86_64") then
    cd centos4.0_x86_64
else 
    cd $PLATFORM
endif

if (( "$RELEASE_TYPE" != "dev") && ( "$RELEASE_TYPE" != "stable-pub")) then
  echo "ERROR: release_type is either dev or stable-pub"
  exit 1
endif

if ("$RELEASE_TYPE" == "stable-pub") setenv RELEASE_TYPE stable3-pub

if ( "$PLATFORM" == "tiger") then
  if (-e /Users/Shared/tmp/$RELEASE_TYPE) \
    rm -Rf /Users/Shared/tmp/$RELEASE_TYPE
  cp -R $RELEASE_TYPE /Users/Shared/tmp
  cd /Users/Shared/tmp
endif

if (-e freesurfer) rm freesurfer
set cmd=(ln -s $RELEASE_TYPE freesurfer)
echo $cmd
$cmd

setenv OSTYPE `uname -s`
if ("$OSTYPE" == "linux") setenv OSTYPE Linux
if ("$OSTYPE" == "Linux") setenv OSTYPE Linux
if ("$OSTYPE" == "darwin") setenv OSTYPE Darwin
if ("$OSTYPE" == "Darwin") setenv OSTYPE Darwin
setenv DATE `date +%Y%m%d`
setenv FILENAME freesurfer-${OSTYPE}-${PLATFORM}-${RELEASE_TYPE}${DATE}-full
setenv  TARNAME ${FILENAME}.tar
echo creating $TARNAME...
tar -X ${SPACE_FREESURFER}/build/scripts/exclude_from_targz -hcvf $TARNAME freesurfer
echo gzipping $TARNAME...
gzip $TARNAME
mv $TARNAME.gz ${SPACE_FREESURFER}/build/pub-releases
chmod g+w ${SPACE_FREESURFER}/build/pub-releases/$TARNAME.gz
if (-e /Users/Shared/tmp/$RELEASE_TYPE) rm -Rf /Users/Shared/tmp/$RELEASE_TYPE

echo md5sum $TARNAME...
md5sum ${SPACE_FREESURFER}/build/pub-releases/$TARNAME.gz

# for the Mac, create_dmg creates the .dmg file from the .tar.gz file

# also, the stable release should be renamed with a version number

