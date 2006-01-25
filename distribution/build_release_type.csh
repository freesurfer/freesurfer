#!/bin/tcsh -f

unsetenv echo
if ($?SET_ECHO_1) set echo=1

set VERSION='$Id: build_release_type.csh,v 1.13 2006/01/25 20:50:55 nicks Exp $'
set HOSTNAME=`hostname -s`
setenv OSTYPE `uname -s`
if ("$OSTYPE" == "linux") setenv OSTYPE Linux
if ("$OSTYPE" == "Linux") setenv OSTYPE Linux
if ("$OSTYPE" == "darwin") setenv OSTYPE Darwin
if ("$OSTYPE" == "Darwin") setenv OSTYPE Darwin

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

# this QTDIR path is also used in the configure command further below
setenv QTDIR /usr/pubsw/packages/qt/current
setenv GLUT_DYLIB_DIR ""
# on Mac OS X Tiger, glut is not automatically in lib path.
# also, need /sw/bin to get latex and dvips
if ("$OSTYPE" == "Darwin") then
  set GLUT_DYLIB_DIR=/usr/pubsw/packages/tiffjpegglut/lib
  setenv PATH "/sw/bin":"$PATH"
  rehash
endif
setenv LD_LIBRARY_PATH "${QTDIR}/lib":"${GLUT_DYLIB_DIR}"
setenv DYLD_LIBRARY_PATH "${QTDIR}/lib":"${GLUT_DYLIB_DIR}"

# Output file
######################################################################
set MAIL_LIST=(kteich@nmr.mgh.harvard.edu nicks@nmr.mgh.harvard.edu)
set FAILURE_MAIL_LIST=(fsdev@nmr.mgh.harvard.edu)
set OUTPUTF=$LOG_DIR/build_log-$RELEASE_TYPE-$HOSTNAME.txt
echo "$HOSTNAME $RELEASE_TYPE build" >& $OUTPUTF
chmod g+w $OUTPUTF
set BEGIN_TIME=`date`
echo $BEGIN_TIME >>& $OUTPUTF

set TIME_STAMP=`date +%Y%m%d`
set PLATFORM=`cat /usr/local/freesurfer/PLATFORM`
set OS=`uname -s`

# Sanity checks
######################################################################
if(! -d $SCRIPT_DIR) then 
  echo "$SCRIPT_DIR doesn't exist" >>& $OUTPUTF
  mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - sanity" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif
if(! -d $DEV_DIR) then 
  echo "$DEV_DIR doesn't exist" >>& $OUTPUTF
  mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - sanity" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif
if(! -d $DEST_DIR) then 
  echo "$DEST_DIR doesn't exist" >>& $OUTPUTF
  mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - sanity" $FAILURE_MAIL_LIST < $OUTPUTF
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
echo "DEST_DIR $DEST_DIR" >>& $OUTPUTF
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
echo "CMD: cd $DEV_DIR" >>& $OUTPUTF
cd ${DEV_DIR} >>& $OUTPUTF
echo "CMD: cvs update -d \>\& $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
cvs update -d >& $LOG_DIR/update-output-$HOSTNAME
chmod g+w $LOG_DIR/update-output-$HOSTNAME

echo "CMD: grep -e "Permission denied" $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
grep -e "Permission denied"  $LOG_DIR/update-output-$HOSTNAME >& /dev/null
if ($status == 0) then
  echo "cvs update: Permission denied" >>& $OUTPUTF
  mail -s "$HOSTNAME $RELEASE_TYPE build FAILED - cvs update permission denied" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e ^\[UP\]\  $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
grep -e ^\[UP\]\  $LOG_DIR/update-output-$HOSTNAME >& /dev/null
if ($status != 0 && ! -e ${FAILED_FILE} ) then
  echo "Nothing changed in repository, SKIPPED building" >>& $OUTPUTF
  mail -s "$HOSTNAME $RELEASE_TYPE build skipped - no cvs changes" $MAIL_LIST < $OUTPUTF
  exit 0
endif

# assume failure (file removed only after successful build)
touch ${FAILED_FILE}
chmod g+w ${FAILED_FILE}

echo "CMD: cat $LOG_DIR/update-output-$HOSTNAME \>\>\& $OUTPUTF" >>& $OUTPUTF
cat $LOG_DIR/update-output-$HOSTNAME >>& $OUTPUTF
echo "CMD: rm -f $LOG_DIR/update-output-$HOSTNAME" >>& $OUTPUTF
rm -f $LOG_DIR/update-output-$HOSTNAME

echo "##########################################################" >>& $OUTPUTF
echo "Freshening Makefiles" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: make distclean" >>& $OUTPUTF
if (-e Makefile) make distclean >>& $OUTPUTF
echo "CMD: rm -rf autom4te.cache" >>& $OUTPUTF
if (-e autom4te.cache) rm -rf autom4te.cache >>& $OUTPUTF
echo "CMD: libtoolize --force" >>& $OUTPUTF
if ( "`uname -s`" == "Linux") libtoolize --force >>& $OUTPUTF
if ( "`uname -s`" == "Darwin") glibtoolize --force >>& $OUTPUTF
echo "CMD: autoreconf --force" >>& $OUTPUTF
autoreconf --force >>& $OUTPUTF
echo "CMD: aclocal" >>& $OUTPUTF
aclocal >>& $OUTPUTF
echo "CMD: autoconf" >>& $OUTPUTF
autoconf >>& $OUTPUTF
echo "CMD: automake" >>& $OUTPUTF
automake >>& $OUTPUTF
echo "CMD: ./configure..." >>& $OUTPUTF
# notice that the configure command sets 'bindir' to /bin-new, overriding
# the default /bin.  later, after make install, bin-new is moved to /bin.
# this is to minimize disruption of machines running recon-all.
./configure \
--with-mni-dir=/usr/pubsw/packages/mni/current \
--with-gsl-dir=/usr/pubsw/packages/gsl/current \
--with-tcl-dir=/usr/pubsw/packages/tcltktixblt/current \
--with-qt-dir=${QTDIR} \
--prefix=${DEST_DIR} \
--bindir=${DEST_DIR}/bin-new \
--enable-nmr-install \
`cat ${BUILD_DIR}/configure_options.txt` >>& $OUTPUTF
if ($status != 0) then
  echo "########################################################" >>& $OUTPUTF
  echo "config.log" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  cat ${DEV_DIR}/config.log >>& $OUTPUTF
  mail -s "$HOSTNAME $RELEASE_TYPE build FAILED after configure" $FAILURE_MAIL_LIST < $OUTPUTF
  touch ${FAILED_FILE}
  chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${DEV_DIR}/*.*" >>& $OUTPUTF
  chmod -R g+rw ${DEV_DIR}/*.* >>& $OUTPUTF
  exit 1
endif

echo "##########################################################" >>& $OUTPUTF
echo "Building dev" >>& $OUTPUTF
echo "" >>& $OUTPUTF
# make
echo "CMD: make" >>& $OUTPUTF
make >>& $OUTPUTF
if ($status != 0) then
  # note: /usr/local/freesurfer/dev/bin/ dirs have not 
  # been modified (bin/ gets written after make install)
  mail -s "$HOSTNAME $RELEASE_TYPE build (make) FAILED" $FAILURE_MAIL_LIST < $OUTPUTF
  touch ${FAILED_FILE}
  chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${DEV_DIR}/*.*" >>& $OUTPUTF
  chmod -R g+rw ${DEV_DIR}/*.* >>& $OUTPUTF
  exit 1  
endif

# make install
# (recall that configure sets $bindir to bin-new/ instead of /bin, 
# to minimize disruption of machines using contents of /bin)
echo "CMD: rm -Rf ${DEST_DIR}/bin-new" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-new) rm -rf ${DEST_DIR}/bin-new >>& $OUTPUTF
echo "CMD: make install" >>& $OUTPUTF
make install >>& $OUTPUTF
if ($status != 0) then
  mail -s "$HOSTNAME $RELEASE_TYPE build (make install) FAILED" $FAILURE_MAIL_LIST < $OUTPUTF
  touch ${FAILED_FILE}
  chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${DEV_DIR}/*.*" >>& $OUTPUTF
  chmod -R g+rw ${DEV_DIR}/*.* >>& $OUTPUTF
  exit 1  
endif

# On the Mac, the Qt-based packages have special targets for installing
# files necessary to run as Mac Apps
if ("$OSTYPE" == "Darwin") then
  echo "CMD: make scuba_bundle" >>& $OUTPUTF
  cd scuba
  make scuba_bundle >>& $OUTPUTF
  echo "CMD: make qdec_bundle" >>& $OUTPUTF
  cd ../qdec
  make qdec_bundle >>& $OUTPUTF
  echo "CMD: make plotter_bundle" >>& $OUTPUTF
  cd ../plotter
  make plotter_bundle >>& $OUTPUTF
  cd ..
endif

# Shift bin/ to bin-old/, and bin-old/ to bin-old-old/ to keep around old versions.
# Move bin/ to bin-old/ instead of copy, to avoid core dumps if some script
# is using a binary in bin/.
# Move newly created bin-new/ to bin/.
# This series of mv's minimizes the time window where the /bin directory
# would appear empty to a machine trying to reference its contents in recon-all
echo "CMD: rm -rf ${DEST_DIR}/bin-old-old" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-old-old) rm -rf ${DEST_DIR}/bin-old-old >>& $OUTPUTF
echo "CMD: mv ${DEST_DIR}/bin-old ${DEST_DIR}/bin-old-old" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-old) mv ${DEST_DIR}/bin-old ${DEST_DIR}/bin-old-old >>& $OUTPUTF
echo "CMD: mv ${DEST_DIR}/bin ${DEST_DIR}/bin-old" >>& $OUTPUTF
mv ${DEST_DIR}/bin ${DEST_DIR}/bin-old >>& $OUTPUTF
echo "CMD: mv ${DEST_DIR}/bin-new ${DEST_DIR}/bin" >>& $OUTPUTF
mv ${DEST_DIR}/bin-new ${DEST_DIR}/bin >>& $OUTPUTF

echo "##########################################################" >>& $OUTPUTF
echo "Setting permissions" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: chmod -R g+rw ${DEST_DIR}/*.*" >>& $OUTPUTF
chmod -R g+rw ${DEST_DIR}/*.* >>& $OUTPUTF
echo "CMD: chmod -R g+rw ${DEV_DIR}/*.*" >>& $OUTPUTF
chmod -R g+rw ${DEV_DIR}/*.* >>& $OUTPUTF
echo "CMD: chmod -R g+rw ${LOG_DIR}/*.*" >>& $OUTPUTF
chmod -R g+rw ${LOG_DIR}/*.* >>& $OUTPUTF

if ($?PUB_DEST_DIR) then
  echo "########################################################" >>& $OUTPUTF
  echo "Building public stable" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  # make
  echo "CMD: make release prefix=$PUB_DEST_DIR" >>& $OUTPUTF
  make release prefix=${PUB_DEST_DIR} >>& $OUTPUTF
  if ($status != 0) then
    mail -s "$HOSTNAME $RELEASE_TYPE release build (make) FAILED" $FAILURE_MAIL_LIST < $OUTPUTF
    touch ${FAILED_FILE}
    chmod g+w ${FAILED_FILE}
    # set group write bit on files changed by make tools:
    echo "CMD: chmod -R g+rw ${PUB_DEST_DIR}/*.*" >>& $OUTPUTF
    chmod -R g+rw ${PUB_DEST_DIR}/*.* >>& $OUTPUTF
    exit 1  
  endif
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${PUB_DEST_DIR}/*.*" >>& $OUTPUTF
  chmod -R g+rw ${PUB_DEST_DIR}/*.* >>& $OUTPUTF
endif

# Success, so remove fail indicator:
rm -rf ${FAILED_FILE}

done:

echo "##########################################################" >>& $OUTPUTF
echo "Done." >>& $OUTPUTF
set END_TIME=`date`
echo $END_TIME >>& $OUTPUTF

# Finish up
######################################################################

# Move log file to stamped version.
chmod g+w $OUTPUTF
mv $OUTPUTF ${LOG_DIR}/build_log-$RELEASE_TYPE-$HOSTNAME-$TIME_STAMP.txt
gzip -f ${LOG_DIR}/build_log-$RELEASE_TYPE-$HOSTNAME-$TIME_STAMP.txt

# Send email.
echo "Begin ${BEGIN_TIME}, end ${END_TIME}" >& $LOG_DIR/message-$HOSTNAME.txt
mail -s "$HOSTNAME $RELEASE_TYPE build is wicked awesome." $MAIL_LIST < $LOG_DIR/message-$HOSTNAME.txt
rm $LOG_DIR/message-$HOSTNAME.txt

# Soon, processor-specific builds:
# make install ${x8664BUILDCXXFLAGS} bindir=\\\$\\\{prefix\\\}/bin/x86-64
# make release prefix={$STABLEDEST} bindir=\\\$\\\{prefix\\\}/bin/x86-64
