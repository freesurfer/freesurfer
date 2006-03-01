#!/bin/tcsh -f

set VERSION='$Id: build_release_type.csh,v 1.36 2006/03/01 00:51:42 nicks Exp $'
unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

set HOSTNAME=`hostname -s`
setenv OSTYPE `uname -s`
if ("$OSTYPE" == "linux") setenv OSTYPE Linux
if ("$OSTYPE" == "Linux") setenv OSTYPE Linux
if ("$OSTYPE" == "darwin") setenv OSTYPE Darwin
if ("$OSTYPE" == "Darwin") setenv OSTYPE Darwin
set OS=${OSTYPE}
set PLATFORM=`cat /usr/local/freesurfer/PLATFORM`

# Set up directories.
######################################################################
#
setenv BUILD_DIR /space/freesurfer/build/$HOSTNAME

if ("$1" == "dev") then
  set RELEASE_TYPE=dev
  set DEV_DIR=${BUILD_DIR}/trunk/dev
  set DEST_DIR=/usr/local/freesurfer/dev
else if ("$1" == "stable-pub") then
  set RELEASE_TYPE=stable
  set DEV_DIR=${BUILD_DIR}/stable/dev
  # notice that the destination is the 'stable3' directory
  set DEST_DIR=/usr/local/freesurfer/stable3
  set PUB_DEST_DIR=/usr/local/freesurfer/stable3-pub
else
  echo "ERROR: release_type must be either dev or stable-pub"
  echo ""
  echo "Examples: "
  echo "  build_release_type dev"
  echo "  build_release_type stable-pub"
  exit 1
endif
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

# Output files (OUTPUTF and CVSUPDATEF)
######################################################################
#
set MAIL_LIST=(kteich@nmr.mgh.harvard.edu nicks@nmr.mgh.harvard.edu)
#set FAILURE_MAIL_LIST=(fsdev@nmr.mgh.harvard.edu)
set FAILURE_MAIL_LIST=(nicks@nmr.mgh.harvard.edu)
set FAILED_FILE=${BUILD_DIR}/${RELEASE_TYPE}-build-FAILED
set OUTPUTF=${LOG_DIR}/build_log-${RELEASE_TYPE}-${HOSTNAME}.txt
set CVSUPDATEF=${LOG_DIR}/update-output-${RELEASE_TYPE}-${HOSTNAME}.txt
echo "$HOSTNAME $RELEASE_TYPE build" >& $OUTPUTF
chmod g+w $OUTPUTF
set BEGIN_TIME=`date`
echo $BEGIN_TIME >>& $OUTPUTF
set TIME_STAMP=`date +%Y%m%d`

# Sanity checks
######################################################################
#
if(! -d $SCRIPT_DIR) then 
  echo "$SCRIPT_DIR doesn't exist" >>& $OUTPUTF
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - sanity"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif
if(! -d $DEV_DIR) then 
  echo "$DEV_DIR doesn't exist" >>& $OUTPUTF
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - sanity"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif
if(! -d $DEST_DIR) then 
  echo "$DEST_DIR doesn't exist" >>& $OUTPUTF
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - sanity"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
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
#
# Go to dev directory, update code, and check the result. If there are
# lines starting with "U " or "P " then we had some changes, so go
# through with the build. If not, quit now. But don't quit if the file
# FAILED exists, because that means that the last build failed.
# Also check for 'Permission denied" and "File is in the way" errors.
# Also check for modified files, which is bad, as this checkout is not
# supposed to be used for development, and it means the real file (the
# one in CVS) will not be used.  Also check for removed files, added
# files, and files with conflicts, all these being a big no-no.
echo "##########################################################" >>& $OUTPUTF
echo "Updating $DEV_DIR" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: cd $DEV_DIR" >>& $OUTPUTF
cd ${DEV_DIR} >>& $OUTPUTF
echo "CMD: cvs update -P -d \>\& $CVSUPDATEF" >>& $OUTPUTF
cvs update -P -d >& $CVSUPDATEF
chmod g+w $CVSUPDATEF

echo "CMD: grep -e "Permission denied" $CVSUPDATEF" >>& $OUTPUTF
grep -e "Permission denied" $CVSUPDATEF >& /dev/null
if ($status == 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - cvs update permission denied"
  echo "$msg" >>& $OUTPUTF
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e "cvs update: move away" $CVSUPDATEF" >>& $OUTPUTF
grep -e "cvs update: move away" $CVSUPDATEF >& /dev/null
if ($status == 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - cvs update: file in the way"
  echo "$msg" >>& $OUTPUTF
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e ^\[M\]\  $CVSUPDATEF" >>& $OUTPUTF
grep -e ^\[M\]\   $CVSUPDATEF >& /dev/null
if ($status == 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - cvs update: file modified!"
  echo "$msg" >>& $OUTPUTF
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e ^\[C\]\  $CVSUPDATEF" >>& $OUTPUTF
grep -e ^\[C\]\   $CVSUPDATEF >& /dev/null
if ($status == 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - cvs update: file conflict!"
  echo "$msg" >>& $OUTPUTF
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e ^\[R\]\  $CVSUPDATEF" >>& $OUTPUTF
grep -e ^\[R\]\   $CVSUPDATEF >& /dev/null
if ($status == 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - cvs update: file removed!"
  echo "$msg" >>& $OUTPUTF
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e ^\[A\]\  $CVSUPDATEF" >>& $OUTPUTF
grep -e ^\[A\]\   $CVSUPDATEF >& /dev/null
if ($status == 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - cvs update: file added!"
  echo "$msg" >>& $OUTPUTF
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif

echo "CMD: grep -e ^\[UP\]\  $CVSUPDATEF" >>& $OUTPUTF
grep -e ^\[UP\]\   $CVSUPDATEF >& /dev/null
if ($status != 0 && ! -e ${FAILED_FILE} ) then
  echo "Nothing changed in repository, SKIPPED building" >>& $OUTPUTF
  set msg="$HOSTNAME $RELEASE_TYPE build skipped - no cvs changes"
  mail -s "$msg" $MAIL_LIST < $OUTPUTF
  echo "CMD: cat $CVSUPDATEF \>\>\& $OUTPUTF" >>& $OUTPUTF
  cat $CVSUPDATEF >>& $OUTPUTF
  echo "CMD: rm -f $CVSUPDATEF" >>& $OUTPUTF
  rm -f $CVSUPDATEF
  exit 0
endif

# assume failure (file removed only after successful build)
touch ${FAILED_FILE}
chmod g+w ${FAILED_FILE}

echo "CMD: cat $CVSUPDATEF \>\>\& $OUTPUTF" >>& $OUTPUTF
cat $CVSUPDATEF >>& $OUTPUTF
echo "CMD: rm -f $CVSUPDATEF" >>& $OUTPUTF
rm -f $CVSUPDATEF

#
# CVS update is now complete, so now, make distclean, and re-configure
#
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
--with-tixwish=/usr/pubsw/packages/tcltktixblt/current/bin/tixwish8.1.8.4 \
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
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED after configure"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  touch ${FAILED_FILE}
  chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
  chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
  chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
  chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
  exit 1
endif

#
# make
#
echo "##########################################################" >>& $OUTPUTF
echo "Making $DEV_DIR" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: make" >>& $OUTPUTF
make >>& $OUTPUTF
if ($status != 0) then
  # note: /usr/local/freesurfer/dev/bin/ dirs have not 
  # been modified (bin/ gets written after make install)
  set msg="$HOSTNAME $RELEASE_TYPE build (make) FAILED"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  touch ${FAILED_FILE}
  chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
  chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
  chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
  chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
  exit 1  
endif

#
# make install
#
# (recall that configure sets $bindir to bin-new/ instead of /bin, 
# to minimize disruption of machines using contents of /bin)
echo "CMD: rm -Rf ${DEST_DIR}/bin-new" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-new) rm -rf ${DEST_DIR}/bin-new >>& $OUTPUTF
echo "CMD: make install" >>& $OUTPUTF
make install >>& $OUTPUTF
if ($status != 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build (make install) FAILED"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  touch ${FAILED_FILE}
  chmod g+w ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
  chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
  chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
  chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
  exit 1  
endif
# strip symbols from binaries, greatly reducing their size
#strip ${DEST_DIR}/bin-new/* >& /dev/null

# Shift bin/ to bin-old/, and bin-old/ to bin-old-old/ to keep old versions.
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

#
# make install is now complete, and /bin dir is now setup with new code
#

echo "##########################################################" >>& $OUTPUTF
echo "Setting permissions" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: chmod -R g+rw ${DEST_DIR}" >>& $OUTPUTF
chmod -R g+rw ${DEST_DIR} >>& $OUTPUTF
echo "CMD: chmod -R g+rw ${DEV_DIR}" >>& $OUTPUTF
chmod -R g+rw ${DEV_DIR} >>& $OUTPUTF
chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
echo "CMD: chmod -R g+rw ${LOG_DIR}" >>& $OUTPUTF
chmod -R g+rw ${LOG_DIR} >>& $OUTPUTF

#
# If building the stable release, then do the special stuff necessary
# for the public version of it.
#
if ($?PUB_DEST_DIR) then
  echo "########################################################" >>& $OUTPUTF
  echo "Building public stable" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  echo "CMD: make release prefix=$PUB_DEST_DIR" >>& $OUTPUTF
  make release prefix=${PUB_DEST_DIR} >>& $OUTPUTF
  if ($status != 0) then
    set msg="$HOSTNAME $RELEASE_TYPE release build (make) FAILED"
    mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
    touch ${FAILED_FILE}
    chmod g+w ${FAILED_FILE}
    # set group write bit on files changed by make tools:
    echo "CMD: chmod -R g+rw ${PUB_DEST_DIR}" >>& $OUTPUTF
    chmod -R g+rw ${PUB_DEST_DIR} >>& $OUTPUTF
    exit 1  
  endif
  mv ${DEST_DIR}/bin-new ${PUB_DEST_DIR}/bin >>& $OUTPUTF
  # strip symbols from binaries, greatly reducing their size
  strip ${PUB_DEST_DIR}/bin/* >& /dev/null
  # set group write bit on files changed by make tools:
  echo "CMD: chmod -R g+rw ${PUB_DEST_DIR}" >>& $OUTPUTF
  chmod -R g+rw ${PUB_DEST_DIR} >>& $OUTPUTF
endif

#
# ensure that the symlinks to the necessary packages are in place
#
symlinks:

set DEST_DIR_LIST=()
if ($?DEST_DIR) then
  set DEST_DIR_LIST=($DEST_DIR_LIST $DEST_DIR)
endif
if ($?PUB_DEST_DIR) then
  set DEST_DIR_LIST=($DEST_DIR_LIST $PUB_DEST_DIR)
endif
foreach destdir ($DEST_DIR_LIST)
  if (! -d $destdir/mni ) then
    set cmd=(ln -s /usr/pubsw/packages/mni/current $destdir/mni)
    echo "$cmd" >>& $OUTPUTF
    $cmd
  endif
  if (! -d $destdir/fsl ) then
    set cmd=(ln -s /usr/pubsw/packages/fsl/current $destdir/fsl)
    echo "$cmd" >>& $OUTPUTF
    $cmd
  endif
  if (! -d $destdir/lib/tcltktixblt ) then
    set cmd=(ln -s /usr/pubsw/packages/tcltktixblt/current $destdir/lib/tcltktixblt)
    echo "$cmd" >>& $OUTPUTF
    $cmd
  endif
  if (! -d $destdir/lib/gsl ) then
    set cmd=(ln -s /usr/pubsw/packages/gsl/current $destdir/lib/gsl)
    echo "$cmd" >>& $OUTPUTF
    $cmd
  endif
  if (! -d $destdir/lib/qt ) then
    set cmd=(ln -s /usr/pubsw/packages/qt/current $destdir/lib/qt)
    echo "$cmd" >>& $OUTPUTF
    $cmd
  endif
  if ("$OSTYPE" == "Darwin") then
    if (! -d $destdir/lib/misc ) then
      set cmd=(ln -s /usr/pubsw/packages/tiffjpegglut $destdir/lib/misc)
      echo "$cmd" >>& $OUTPUTF
      $cmd
    endif
  endif
  # also sample subject:
  if (! -d $destdir/subjects/bert ) then
    set cmd=(ln -s /space/freesurfer/subjects/bert $destdir/subjects/bert)
    echo "$cmd" >>& $OUTPUTF
    $cmd
  endif
end


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
set msg="$HOSTNAME $RELEASE_TYPE build is wicked awesome."
mail -s "$msg" $MAIL_LIST < $LOG_DIR/message-$HOSTNAME.txt
rm $LOG_DIR/message-$HOSTNAME.txt
