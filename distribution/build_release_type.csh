#!/bin/tcsh -f

set ECHO=
#set echo=1

set HOSTNAME=`hostname -s`

# Set up directories.
######################################################################
if(! $?BUILD_DIR) then
  setenv BUILD_DIR /space/freesurfer/build/$HOSTNAME
endif

if ("$1" == "dev") then
  set RELEASE_TYPE=dev
  set DEV_DIR=${BUILD_DIR}/trunk/dev
  set DEST_DIR=/usr/local/freesurfer/dev
else if ("$1" == "stable") then
  set RELEASE_TYPE=stable
  set DEV_DIR=${BUILD_DIR}/stable/dev
  set DEST_DIR=/usr/local/freesurfer/beta
  set PUB_DEST_DIR=/usr/local/freesurfer/pub
else
  echo "ERROR: release_type must be either dev or stable"
  exit 1
endif

set FAILED_FILE=${BUILD_DIR}/${RELEASE_TYPE}-build-FAILED
set SCRIPT_DIR=/space/freesurfer/build/scripts
set LOG_DIR=/space/freesurfer/build/logs

# QTDIR should match that declared in each $HOSTNAME/configure_options.txt file
setenv QTDIR /usr/pubsw/packages/qt/current
setenv GLUT_DYLIB_DIR ""
# on Mac OS X Tiger, glut is not automatically in lib path.
# also, need /sw/bin to get latex and dvips
if ("$HOSTNAME" == "storm") then
    set GLUT_DYLIB_DIR=/usr/pubsw/packages/tiffjpegglut/lib
    setenv PATH "/sw/bin":"$PATH"
    rehash
endif
setenv LD_LIBRARY_PATH "${QTDIR}/lib":"${GLUT_DYLIB_DIR}"
setenv DYLD_LIBRARY_PATH "${QTDIR}/lib":"${GLUT_DYLIB_DIR}"

# Output file
######################################################################
set MAIL_LIST=(kteich@nmr.mgh.harvard.edu nicks@nmr.mgh.harvard.edu)
set OUTPUTF=$LOG_DIR/build_log-$RELEASE_TYPE-$HOSTNAME.txt
echo "$HOSTNAME $RELEASE_TYPE build" >& $OUTPUTF
$ECHO chmod g+w $OUTPUTF
set BEGIN_TIME=`date`
echo $BEGIN_TIME >>& $OUTPUTF

set TIME_STAMP=`date +%Y%m%d`
set PLATFORM=`cat /usr/local/freesurfer/PLATFORM`
set OS=`uname -s`

# Sanity checks
######################################################################
if(! -d $SCRIPT_DIR) then 
  echo "$SCRIPT_DIR doesn't exist" >>& $OUTPUTF
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - sanity" $MAIL_LIST < $OUTPUTF
  exit 1  
endif
if(! -d $DEV_DIR) then 
  echo "$DEV_DIR doesn't exist" >>& $OUTPUTF
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - sanity" $MAIL_LIST < $OUTPUTF
  exit 1  
endif

# processor-specific build options.
######################################################################
#set P3CXXFLAGS="CXXFLAGS=-march=pentium3"
#set P4CXXFLAGS="CXXFLAGS=-march=pentium4-64"
#set x8664CXXFLAGS="CXXFLAGS=-march=x86-64"

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
echo "DEV_DIR $DEV_DIR" >>& $OUTPUTF
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
$ECHO echo "CMD: cvs update -d \>\& $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
$ECHO cvs update -d >& $LOG_DIR/update-output-$HOSTNAME

$ECHO echo "CMD: grep -e "Permission denied" $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
$ECHO grep -e "Permission denied"  $LOG_DIR/update-output-$HOSTNAME >& /dev/null
if ($status == 0) then
  echo "cvs update: Permission denied" >>& $OUTPUTF
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - cvs update permission denied" $MAIL_LIST < $OUTPUTF
  exit 1  
endif

$ECHO echo "CMD: grep -e ^\[UP\]\  $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
$ECHO grep -e ^\[UP\]\  $LOG_DIR/update-output-$HOSTNAME >& /dev/null
if ($status != 0 && ! -e ${FAILED_FILE} ) then
  echo "Nothing changed in repository, SKIPPED building" >>& $OUTPUTF
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build skipped - no cvs changes" $MAIL_LIST < $OUTPUTF
  exit 0
endif

# assume failure (file removed only after successful build)
$ECHO touch ${FAILED_FILE}
$ECHO chmod g+w ${FAILED_FILE}

$ECHO echo "CMD: cat $LOG_DIR/update-output-$HOSTNAME \>\>\& $OUTPUTF" >>& $OUTPUTF
$ECHO cat $LOG_DIR/update-output-$HOSTNAME >>& $OUTPUTF
$ECHO echo "CMD: rm -f $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
$ECHO rm -f $LOG_DIR/update-output-$HOSTNAME

echo "##########################################################" >>& $OUTPUTF
echo "Freshening Makefiles" >>& $OUTPUTF
echo "" >>& $OUTPUTF
$ECHO echo "CMD: make distclean" >>& $OUTPUTF
if (-e Makefile) make distclean >>& $OUTPUTF
$ECHO echo "CMD: rm -rf autom4te.cache" >>& $OUTPUTF
if (-e autom4te.cache) rm -rf autom4te.cache >>& $OUTPUTF
$ECHO echo "CMD: libtoolize --force" >>& $OUTPUTF
if ( "`uname -s`" == "Linux") libtoolize --force >>& $OUTPUTF
if ( "`uname -s`" == "Darwin") glibtoolize --force >>& $OUTPUTF
$ECHO echo "CMD: autoreconf --force" >>& $OUTPUTF
$ECHO autoreconf --force >>& $OUTPUTF
$ECHO echo "CMD: aclocal" >>& $OUTPUTF
$ECHO aclocal >>& $OUTPUTF
$ECHO echo "CMD: autoconf" >>& $OUTPUTF
$ECHO autoconf >>& $OUTPUTF
$ECHO echo "CMD: automake" >>& $OUTPUTF
$ECHO automake >>& $OUTPUTF
$ECHO echo "CMD: ./configure `cat ${BUILD_DIR}/configure_options.txt` --prefix=${DEST_DIR} --enable-nmr-install" >>& $OUTPUTF
$ECHO ./configure `cat ${BUILD_DIR}/configure_options.txt` --prefix=${DEST_DIR} --enable-nmr-install >>& $OUTPUTF
if ($status != 0) then
  echo "########################################################" >>& $OUTPUTF
  echo "config.log" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  $ECHO cat ${DEV_DIR}/config.log >>& $OUTPUTF
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build FAILED after configure" $MAIL_LIST < $OUTPUTF
  $ECHO touch ${FAILED_FILE}
  $ECHO chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  $ECHO echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
  $ECHO chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
  exit 1  
endif

echo "##########################################################" >>& $OUTPUTF
echo "Building dev" >>& $OUTPUTF
echo "" >>& $OUTPUTF
# make
$ECHO echo "CMD: make -j 4" >>& $OUTPUTF
$ECHO make -j 4 >>& $OUTPUTF
if ($status != 0) then
  # note: /usr/local/freesurfer/dev/bin/ dirs have not 
  # been modified (make install does that)
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build (make) FAILED" $MAIL_LIST < $OUTPUTF
  $ECHO touch ${FAILED_FILE}
  $ECHO chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  $ECHO echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
  $ECHO chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
  exit 1  
endif

# Shift bin to bin-old and bin-old to bin-old-old to keep around old versions.
$ECHO echo "CMD: rm -rf ${DEST_DIR}/bin-old-old" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-old-old) rm -rf ${DEST_DIR}/bin-old-old >>& $OUTPUTF
$ECHO echo "CMD: mv ${DEST_DIR}/bin-old ${DEST_DIR}/bin-old-old" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-old) mv ${DEST_DIR}/bin-old ${DEST_DIR}/bin-old-old >>& $OUTPUTF
$ECHO echo "CMD: mv ${DEST_DIR}/bin ${DEST_DIR}/bin-old" >>& $OUTPUTF
$ECHO mv ${DEST_DIR}/bin ${DEST_DIR}/bin-old >>& $OUTPUTF

# make install (-j 4 flag launches up to 4 jobs simultaneously, using multi-processors)
$ECHO echo "CMD: make -j 4 install" >>& $OUTPUTF
$ECHO make -j 4 install >>& $OUTPUTF
if ($status != 0) then
  $ECHO mail -s "$HOSTNAME $RELEASE_TYPE build (make install) FAILED" $MAIL_LIST < $OUTPUTF
  $ECHO touch ${FAILED_FILE}
  $ECHO chmod g+w ${FAILED_FILE}
  # restore prior (possibly working) bin/ dirs
  $ECHO echo "CMD: mv ${DEST_DIR}/bin ${DEST_DIR}/bin-failed" >>& $OUTPUTF
  $ECHO mv ${DEST_DIR}/bin ${DEST_DIR}/bin-failed >>& $OUTPUTF
  $ECHO echo "CMD: mv ${DEST_DIR}/bin-old ${DEST_DIR}/bin" >>& $OUTPUTF
  $ECHO mv ${DEST_DIR}/bin-old ${DEST_DIR}/bin >>& $OUTPUTF
  $ECHO echo "CMD: mv ${DEST_DIR}/bin-old-old ${DEST_DIR}/bin-old" >>& $OUTPUTF
  $ECHO mv ${DEST_DIR}/bin-old-old ${DEST_DIR}/bin-old >>& $OUTPUTF
  # set group write bit on files changed by make tools:
  $ECHO echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
  $ECHO chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
  exit 1  
endif

echo "##########################################################" >>& $OUTPUTF
echo "Setting permissions" >>& $OUTPUTF
echo "" >>& $OUTPUTF
$ECHO echo "CMD: chmod -R g+rw ${DEST_DIR}" >>& $OUTPUTF
$ECHO chmod -R g+rw ${DEST_DIR} >>& $OUTPUTF
$ECHO echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
$ECHO chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
$ECHO echo "CMD: chmod -R g+rw ${LOG_DIR}" >>& $OUTPUTF
$ECHO chmod -R g+rw ${LOG_DIR} >>& $OUTPUTF

if ($?PUB_DEST_DIR) then
  echo "########################################################" >>& $OUTPUTF
  echo "Building public stable" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  # make
  $ECHO echo "CMD: make -j 4 release prefix=$PUB_DEST_DIR" >>& $OUTPUTF
  $ECHO make -j 4 release prefix=${PUB_DEST_DIR} >>& $OUTPUTF
  if ($status != 0) then
    $ECHO mail -s "$HOSTNAME $RELEASE_TYPE release build (make) FAILED" $MAIL_LIST < $OUTPUTF
    $ECHO touch ${FAILED_FILE}
    $ECHO chmod g+w ${FAILED_FILE}
    # set group write bit on files changed by make tools:
    $ECHO echo "CMD: chmod -R g+rw ${PUB_DEST_DIR}" >>& $OUTPUTF
    $ECHO chmod -R g+rw ${PUB_DEST_DIR} >>& $OUTPUTF
    exit 1  
  endif
  # set group write bit on files changed by make tools:
  $ECHO echo "CMD: chmod -R g+rw ${PUB_DEST_DIR}" >>& $OUTPUTF
  $ECHO chmod -R g+rw ${PUB_DEST_DIR} >>& $OUTPUTF
endif

# Success, so remove fail indicator:
$ECHO rm -rf ${FAILED_FILE}

done:

echo "##########################################################" >>& $OUTPUTF
echo "Done." >>& $OUTPUTF
set END_TIME=`date`
echo $END_TIME >>& $OUTPUTF

# Finish up
######################################################################

# Move log file to stamped version.
$ECHO mv $OUTPUTF ${LOG_DIR}/build_log-$RELEASE_TYPE-$HOSTNAME-$TIME_STAMP.txt
$ECHO gzip -f ${LOG_DIR}/build_log-$RELEASE_TYPE-$HOSTNAME-$TIME_STAMP.txt

# Send email.
echo "Begin ${BEGIN_TIME}, end ${END_TIME}" >& $LOG_DIR/message-$HOSTNAME.txt
$ECHO mail -s "$HOSTNAME $RELEASE_TYPE build is wicked awesome." $MAIL_LIST < $LOG_DIR/message-$HOSTNAME.txt
$ECHO rm $LOG_DIR/message-$HOSTNAME.txt

# Soon, processor-specific builds:
# $ECHO make install ${x8664BUILDCXXFLAGS} bindir=\\\$\\\{prefix\\\}/bin/x86-64
# $ECHO make install ${x8664BUILDCXXFLAGS} bindir=\\\$\\\{prefix\\\}/bin/x86-64
# $ECHO make release prefix={$STABLEDEST} bindir=\\\$\\\{prefix\\\}/bin/x86-64
