#!/bin/tcsh -ef

unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

if ("$1" == "rh7.3") then
    if ("`uname -n`" != "martinos01" ) then
        echo "must run on machine martinos01"
        exit 1
    endif
else if ("$1" == "rh9") then
    if ("`uname -n`" != "kani" ) then
        echo "must run on machine kani"
        exit 1
    endif
else if ("$1" == "rhel4") then
    if ("`uname -n`" != "triassic" ) then
        echo "must run on machine triassic"
        exit 1
    endif
else if ("$1" == "centos4") then
    if ("`uname -n`" != "fishie" ) then
        echo "must run on machine fishie"
        exit 1
    endif
else if ("$1" == "centos4_x86_64") then
    if ("`uname -n`" != "minerva" ) then
        echo "must run on machine minerva"
        exit 1
    endif
else if ("$1" == "tiger") then
    if ("`uname -n`" != "storm.nmr.mgh.harvard.edu" ) then
        echo "must run on machine storm"
        exit 1
    endif
else
    echo "Usage:"
    echo "$0 <platform> <release_type>"
    echo "where <platform> is rh7.3, rh9, rhel4, centos4, centos4_x86_64, or tiger"
    echo "and <release_type> is either dev, stable or pub"
    exit 1
endif

setenv SPACE_FREESURFER /space/freesurfer
cd $SPACE_FREESURFER
if ( "$1" == "centos4") then
    cd centos4.0
else if ( "$1" == "centos4_x86_64") then
    cd centos4.0_x86_64
else 
    cd $1
endif

if ( "$2" == "dev") then
else if ( "$2" == "stable") then
else if ( "$2" == "pub") then
else 
  echo "ERROR: release_type is either dev, stable or pub"
  exit 1
endif

if ( "$1" == "tiger") then
  if (-e /Users/Shared/tmp/$2) rm -Rf /Users/Shared/tmp/$2
  cp -R $2 /Users/Shared/tmp
  cd /Users/Shared/tmp
endif

if (-e freesurfer) rm freesurfer
ln -s $2 freesurfer

setenv OSTYPE `uname -s`
if ("$OSTYPE" == "linux") setenv OSTYPE Linux
if ("$OSTYPE" == "Linux") setenv OSTYPE Linux
if ("$OSTYPE" == "darwin") setenv OSTYPE Darwin
if ("$OSTYPE" == "Darwin") setenv OSTYPE Darwin
setenv DATE `date +%Y%m%d`
setenv FILENAME freesurfer-${OSTYPE}-$1-${2}${DATE}-full
setenv  TARNAME ${FILENAME}.tar
echo creating $TARNAME...
tar -X ${SPACE_FREESURFER}/build/scripts/exclude_from_targz -hcvf $TARNAME freesurfer
echo gzipping $TARNAME...
gzip $TARNAME
mv $TARNAME.gz ${SPACE_FREESURFER}/build/pub-releases
chmod g+w ${SPACE_FREESURFER}/build/pub-releases/$TARNAME.gz
if (-e /Users/Shared/tmp/$2) rm -Rf /Users/Shared/tmp/$2

echo md5sum $TARNAME...
md5sum ${SPACE_FREESURFER}/build/pub-releases/$TARNAME...

# for the Mac, create_dmg creates the .dmg file from the .tar.gz file

# also, the stable release should be renamed with a version number

