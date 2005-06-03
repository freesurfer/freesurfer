#!/bin/tcsh -f


set ECHO=

# Set up directories.
######################################################################
if(! $?BUILD_DIR) then
  setenv BUILD_DIR /space/birn/50/freesurfer/build/`hostname -s`
endif

setenv QTDIR /space/birn/50/freesurfer/build/`hostname -s`/qt
setenv LD_LIBRARY_PATH ${QTDIR}/lib
setenv DYLD_LIBRARY_PATH ${QTDIR}/lib

set SCRIPT_DIR=/space/birn/50/freesurfer/build/scripts
set LOG_DIR=/space/birn/50/freesurfer/build/logs
set DEV_DIR=${BUILD_DIR}/trunk/dev
set DEV_DEST_DIR=/usr/local/freesurfer/dev

set FAILED_FILE=${BUILD_DIR}/dev-FAILED

# Output file
######################################################################
set MAIL_LIST=(kteich@nmr.mgh.harvard.edu)
set OUTPUTF=/tmp/build_log-dev-`hostname -s`.txt
echo "`hostname -s` build" >& $OUTPUTF
set BEGIN_TIME=`date`
echo $BEGIN_TIME >>& $OUTPUTF

set TIME_STAMP=`date +%Y%m%d`
set PLATFORM=`cat /usr/local/freesurfer/PLATFORM`
set OS=`uname -s`

# Sanity checks
######################################################################
if(! -d $SCRIPT_DIR) then 
  echo "$SCRIPT_DIR doesn't exist" >>& $OUTPUTF
  $ECHO mail -s "`hostname -s` dev build FAILED - sanity" $MAIL_LIST < $OUTPUTF
  exit 1  
endif

if(! -d $DEV_DIR) then 
  echo "$DEV_DIR doesn't exist" >>& $OUTPUTF
  $ECHO mail -s "`hostname -s` dev build FAILED - sanity" $MAIL_LIST < $OUTPUTF
  exit 1  
endif

# Build options.
######################################################################
set P3CXXFLAGS="CXXFLAGS=-march=pentium3"
set P4CXXFLAGS="CXXFLAGS=-march=pentium4-64"
set x8664CXXFLAGS="CXXFLAGS=-march=x86-64"

# Source the source_before_building file if they have it
if( -f ${BUILD_DIR}/source_before_building.csh ) then
  source ${BUILD_DIR}/source_before_building.csh
endif

echo "##########################################################" >>& $OUTPUTF
echo "Settings" >>& $OUTPUTF
echo "BUILD_DIR $BUILD_DIR" >>& $OUTPUTF
echo "QTDIR $QTDIR" >>& $OUTPUTF
echo "LD_LIBRARY_PATH $LD_LIBRARY_PATH" >>& $OUTPUTF
echo "DYLD_LIBRARY_PATH $DYLD_LIBRARY_PATH" >>& $OUTPUTF
echo "SCRIPT_DIR $SCRIPT_DIR" >>& $OUTPUTF
echo "LOG_DIR $LOG_DIR" >>& $OUTPUTF
if( $?CFLAGS ) then 
  echo "CFLAGS $CFLAGS" >>& $OUTPUTF
endif
if( $?CPPFLAGS ) then 
  echo "CPPFLAGS $CPPFLAGS" >>& $OUTPUTF
endif
if( $?CXXFLAGS ) then 
  echo "CXXFLAGS $CXXFLAGS" >>& $OUTPUTF
endif
if( $?LDFLAGS ) then 
  echo "LDFLAGS $LDFLAGS" >>& $OUTPUTF
endif
echo "" >>& $OUTPUTF

# Do the build.
######################################################################

# Go to dev directory, update code, and check the result. If there are
# lines starting with "U " or "P " then we had some changes, so go
# through with the build. If not, quit now. But don't quit if the file
# FAILED exists, because that means that the last build failed.
echo "##########################################################" >>& $OUTPUTF
echo "Updating dev" >>& $OUTPUTF
echo "" >>& $OUTPUTF
$ECHO echo "CMD: cd $DEV_DIR" >>& $OUTPUTF
$ECHO cd $DEV_DIR >>& $OUTPUTF
$ECHO echo "CMD: cvs update -d \>\& /tmp/update-output" >>& $OUTPUTF
$ECHO cvs update -d >& /tmp/update-output

$ECHO echo "CMD: grep -e ^\[UP\]\  /tmp/update-output" >>& $OUTPUTF
$ECHO grep -e ^\[UP\]\  /tmp/update-output >& /dev/null
if ($status != 0 && ! -e ${FAILED_FILE} ) then
  echo "Nothing changed in repository, SKIPPED building" >>& $OUTPUTF
  goto done
endif

$ECHO echo "CMD: cat /tmp/update-output \>\>\& $OUTPUTF" >>& $OUTPUTF
$ECHO cat /tmp/update-output >>& $OUTPUTF
$ECHO echo "CMD: rm -f /tmp/update-output" >>& $OUTPUTF
$ECHO rm -f /tmp/update-output



echo "##########################################################" >>& $OUTPUTF
echo "Freshening Makefiles" >>& $OUTPUTF
echo "" >>& $OUTPUTF
$ECHO echo "CMD: make distclean" >>& $OUTPUTF
$ECHO make distclean >>& $OUTPUTF
$ECHO echo "CMD: rm -rf autom4te.cache" >>& $OUTPUTF
$ECHO rm -rf autom4te.cache >>& $OUTPUTF
$ECHO echo "CMD: aclocal" >>& $OUTPUTF
$ECHO aclocal >>& $OUTPUTF
$ECHO echo "CMD: autoconf" >>& $OUTPUTF
$ECHO autoconf >>& $OUTPUTF
$ECHO echo "CMD: automake" >>& $OUTPUTF
$ECHO automake >>& $OUTPUTF
$ECHO echo "CMD: ./configure `cat ${BUILD_DIR}/configure_options.txt` --prefix=${DEV_DEST_DIR}" >>& $OUTPUTF
$ECHO ./configure `cat ${BUILD_DIR}/configure_options.txt` --prefix=${DEV_DEST_DIR} >>& $OUTPUTF
if ($status != 0) then
  echo "########################################################" >>& $OUTPUTF
  echo "config.log" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  $ECHO cat ${DEV_DIR}/config.log >>& $OUTPUTF
  $ECHO mail -s "`hostname -s` dev build FAILED after dev configure" $MAIL_LIST < $OUTPUTF
  $ECHO touch ${FAILED_FILE}
  exit 1  
endif

echo "##########################################################" >>& $OUTPUTF
echo "Building dev" >>& $OUTPUTF
echo "" >>& $OUTPUTF
$ECHO echo "CMD: rm -rf ${DEV_DEST_DIR}/bin-old" >>& $OUTPUTF
$ECHO rm -rf ${DEV_DEST_DIR}/bin-old >>& $OUTPUTF
$ECHO echo "CMD: cp -r ${DEV_DEST_DIR}/bin ${DEV_DEST_DIR}/bin-old" >>& $OUTPUTF
$ECHO cp -r ${DEV_DEST_DIR}/bin ${DEV_DEST_DIR}/bin-old >>& $OUTPUTF
$ECHO echo "CMD: make install" >>& $OUTPUTF
$ECHO make install >>& $OUTPUTF
if ($status != 0) then
  $ECHO mail -s "`hostname -s` dev build FAILED after building dev" $MAIL_LIST < $OUTPUTF
  $ECHO touch ${FAILED_FILE}
  exit 1  
endif

$ECHO rm -rf ${FAILED_FILE}

done:

echo "##########################################################" >>& $OUTPUTF
echo "Done." >>& $OUTPUTF
set END_TIME=`date`
echo $END_TIME >>& $OUTPUTF

# Finish up
######################################################################

# Move log file to stamped version.
$ECHO mv $OUTPUTF ${LOG_DIR}/build_log-dev-`hostname -s`-$TIME_STAMP.txt
$ECHO gzip -f ${LOG_DIR}/build_log-dev-`hostname -s`-$TIME_STAMP.txt

# Send email.
echo "Begin ${BEGIN_TIME}, end ${END_TIME}" >& /tmp/message.txt
$ECHO mail -s "`hostname -s` dev build is awesome." $MAIL_LIST < /tmp/message.txt
$ECHO rm /tmp/message.txt

# Soon:
# $ECHO make install ${x8664BUILDCXXFLAGS} bindir=\\\$\\\{prefix\\\}/bin/x86-64
# $ECHO make install ${x8664BUILDCXXFLAGS} bindir=\\\$\\\{prefix\\\}/bin/x86-64
# $ECHO make release prefix={$STABLEDEST} bindir=\\\$\\\{prefix\\\}/bin/x86-64
