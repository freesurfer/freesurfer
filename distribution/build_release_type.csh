#!/bin/tcsh -f

set ID='$Id: build_release_type.csh,v 1.87 2007/06/19 15:01:23 nicks Exp $'

unsetenv echo
if ($?SET_ECHO_1) set echo=1

umask 002

# usage:
#  build_release_type dev
#  build_release_type stable
#  build_release_type stable-pub
set RELEASE_TYPE=$1

set STABLE_VER_NUM="v3.0.5"
set STABLE_PUB_VER_NUM="v3.0.5"

set SUCCESS_MAIL_LIST=(nicks kteich)
set FAILURE_MAIL_LIST=(nicks kteich fischl greve dsjen)
#set FAILURE_MAIL_LIST=($SUCCESS_MAIL_LIST)

set HOSTNAME=`hostname -s`
setenv OSTYPE `uname -s`
if ("$OSTYPE" == "linux") setenv OSTYPE Linux
if ("$OSTYPE" == "Linux") setenv OSTYPE Linux
if ("$OSTYPE" == "darwin") setenv OSTYPE Darwin
if ("$OSTYPE" == "Darwin") setenv OSTYPE Darwin
set OS=${OSTYPE}

if ("$OSTYPE" == "Darwin") then
  # Mac OS X chmod and chgrp need -L flag to follow symbolic links
  set change_flags=(-RL)
else
  set change_flags=(-R)
endif

# on minerva, use gcc v4.1
if ("$HOSTNAME" == "minerva") then
#  setenv PATH "/space/minerva/1/users/nicks/pkgs/gcc4.1/install/bin":"$PATH"
endif

#
# Set up directories.
######################################################################
#
setenv SPACE_FS /space/freesurfer
setenv LOCAL_FS /usr/local/freesurfer
# if /space/freesurfer is down, or if there is a need to install
# outside of /usr/local/freesurfer, 
# then the var USE_SPACE_MINERVA can be set
if ($?USE_SPACE_MINERVA) then
  setenv SPACE_FS /space/minerva/1/users/nicks
  setenv LOCAL_FS /space/minerva/1/users/nicks/build/install/${HOSTNAME}
endif

setenv BUILD_DIR      ${SPACE_FS}/build/$HOSTNAME
setenv PLATFORM       "`cat ${LOCAL_FS}/PLATFORM`"
setenv BUILD_PLATFORM "`cat ${BUILD_DIR}/PLATFORM`"

# DEV_DIR is the CVS checkout of either the main trunk and the stable branch.
# DEST_DIR is where the build will be installed.
if ("$RELEASE_TYPE" == "dev") then
  set DEV_DIR=${BUILD_DIR}/trunk/dev
  set DEST_DIR=${LOCAL_FS}/dev
else if ("$RELEASE_TYPE" == "stable") then
  set DEV_DIR=${BUILD_DIR}/stable/dev
  set DEST_DIR=${LOCAL_FS}/stable
else if ("$RELEASE_TYPE" == "stable-pub") then
  set DEV_DIR=${BUILD_DIR}/stable/dev
  set DEST_DIR=${LOCAL_FS}/stable-pub
else
  echo "ERROR: release_type must be either dev, stable or stable-pub"
  echo ""
  echo "Examples: "
  echo "  build_release_type dev"
  echo "  build_release_type stable"
  echo "  build_release_type stable-pub"
  exit 1
endif
set SCRIPT_DIR=${SPACE_FS}/build/scripts
set LOG_DIR=${SPACE_FS}/build/logs

# dev build use latest-and-greatest package libs
# stable build use explicit package versions (for stability)
if (("${RELEASE_TYPE}" == "stable") || ("${RELEASE_TYPE}" == "stable-pub")) then
  set MNIDIR=/usr/pubsw/packages/mni/1.4
  set GSLDIR=/usr/pubsw/packages/gsl/1.6
  set TCLDIR=/usr/pubsw/packages/tcltktixblt/8.4.6
  set TIXWISH=${TCLDIR}/bin/tixwish8.1.8.4
  set VXLDIR=/usr/pubsw/packages/vxl/1.6.0
  set TJGDIR=/usr/pubsw/packages/tiffjpegglut/1.1
  unsetenv QTDIR
  unsetenv FSLDIR
  if (-e /usr/pubsw/packages/fsl/3.2b) then
    setenv FSLDIR /usr/pubsw/packages/fsl/3.2b
  else if (-e /usr/pubsw/packages/fsl/3.2) then
    setenv FSLDIR /usr/pubsw/packages/fsl/3.2
  endif
  unset CPPUNITDIR
else
  # dev build uses most current
  set MNIDIR=/usr/pubsw/packages/mni/current
  set VXLDIR=/usr/pubsw/packages/vxl/current
  set VTKDIR=/usr/pubsw/packages/vtk/current
  set TJGDIR=/usr/pubsw/packages/tiffjpegglut/current
  set TCLDIR=/usr/pubsw/packages/tcltktixblt/current
  set TIXWISH=${TCLDIR}/bin/tixwish8.1.8.4
  setenv FSLDIR /usr/pubsw/packages/fsl/current
  set CPPUNITDIR=/usr/pubsw/packages/cppunit/current
  if ( ! -d ${CPPUNITDIR} ) unset CPPUNITDIR
  # GSL and Qt are no longer used, so they're not defined
  unsetenv QTDIR
  unsetenv GSLDIR
endif

# on Mac OS X Tiger, need /sw/bin (Fink) to get latex and dvips.
if ("$OSTYPE" == "Darwin") then
  setenv PATH "/sw/bin":"$PATH"
  rehash
endif


#
# Output log files (OUTPUTF and CVSUPDATEF)
######################################################################
#
set FAILED_FILE=${BUILD_DIR}/${RELEASE_TYPE}-build-FAILED
set OUTPUTF=${LOG_DIR}/build_log-${RELEASE_TYPE}-${HOSTNAME}.txt
set CVSUPDATEF=${LOG_DIR}/update-output-${RELEASE_TYPE}-${HOSTNAME}.txt
echo "$HOSTNAME $RELEASE_TYPE build" >& $OUTPUTF
chmod g+w $OUTPUTF
set BEGIN_TIME=`date`
echo $BEGIN_TIME >>& $OUTPUTF
set TIME_STAMP=`date +%Y%m%d`

#goto symlinks


#
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
if ("${BUILD_PLATFORM}" != "${PLATFORM}") then
  echo "PLATFORM mismatch!" >>& $OUTPUTF
  echo "${LOCAL_FS}/PLATFORM=${PLATFORM}" >>& $OUTPUTF
  echo "${BUILD_DIR}/PLATFORM=${BUILD_PLATFORM}" >>& $OUTPUTF
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED - sanity"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  exit 1  
endif


# Source the source_before_building file if it exists
if( -f ${BUILD_DIR}/source_before_building.csh ) then
  echo "source ${BUILD_DIR}/source_before_building.csh" >>& $OUTPUTF
  source ${BUILD_DIR}/source_before_building.csh
endif
echo "##########################################################" >>& $OUTPUTF
echo "Settings" >>& $OUTPUTF
echo "PLATFORM $PLATFORM" >>& $OUTPUTF
echo "HOSTNAME $HOSTNAME" >>& $OUTPUTF
echo "BUILD_DIR $BUILD_DIR" >>& $OUTPUTF
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

# in case a new dev dir needed to be created, and the old one cannot
# be deleted because of permissions, then name that old dir 'devold',
# and it will get deleted here:
set DEVOLD=${BUILD_DIR}/trunk/devold
if (-e ${DEVOLD}) then
  echo "CMD: rm -Rf ${DEVOLD}" >>& $OUTPUTF
  rm -rf ${DEVOLD} >>& $OUTPUTF
endif
set DEVOLD=${BUILD_DIR}/stable/devold
if (-e ${DEVOLD}) then
  echo "CMD: rm -Rf ${DEVOLD}" >>& $OUTPUTF
  rm -rf ${DEVOLD} >>& $OUTPUTF
endif


#
# CVS update
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
  mail -s "$msg" $SUCCESS_MAIL_LIST < $OUTPUTF
  echo "CMD: cat $CVSUPDATEF \>\>\& $OUTPUTF" >>& $OUTPUTF
  cat $CVSUPDATEF >>& $OUTPUTF
  echo "CMD: rm -f $CVSUPDATEF" >>& $OUTPUTF
  rm -f $CVSUPDATEF
  exit 0
endif

# assume failure (file removed only after successful build)
rm -f ${FAILED_FILE}
touch ${FAILED_FILE}

echo "CMD: cat $CVSUPDATEF \>\>\& $OUTPUTF" >>& $OUTPUTF
cat $CVSUPDATEF >>& $OUTPUTF
echo "CMD: rm -f $CVSUPDATEF" >>& $OUTPUTF
rm -f $CVSUPDATEF
# CVS update is now complete


#
# make distclean and configure
######################################################################
#
echo "##########################################################" >>& $OUTPUTF
echo "Freshening Makefiles" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: make distclean" >>& $OUTPUTF
if (-e Makefile) make distclean >>& $OUTPUTF
echo "CMD: rm -rf autom4te.cache" >>& $OUTPUTF
if (-e autom4te.cache) rm -rf autom4te.cache >>& $OUTPUTF
echo "CMD: libtoolize --force" >>& $OUTPUTF
if ( "${OSTYPE}" == "Linux") libtoolize --force >>& $OUTPUTF
if ( "${OSTYPE}" == "Darwin") glibtoolize --force >>& $OUTPUTF
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
set ENAB_NMR="--enable-nmr-install"
if ("${RELEASE_TYPE}" == "stable-pub") then
  # public build doesn't get the extra special stuff
  set ENAB_NMR=""
endif
setenv FREESURFER_HOME $DEST_DIR
set cnfgr=(./configure)
set cnfgr=($cnfgr --prefix=${FREESURFER_HOME})
set cnfgr=($cnfgr --bindir=${DEST_DIR}/bin-new)
set cnfgr=($cnfgr $ENAB_NMR)
set cnfgr=($cnfgr `cat ${BUILD_DIR}/configure_options.txt`)
set cnfgr=($cnfgr --with-mni-dir=${MNIDIR})
if ($?GSLDIR) then
    set cnfgr=($cnfgr --with-gsl-dir=${GSLDIR})
endif
set cnfgr=($cnfgr --with-vxl-dir=${VXLDIR})
if ($?VTKDIR) then
    set cnfgr=($cnfgr --with-vtk-dir=${VTKDIR})
endif
set cnfgr=($cnfgr --with-tiffjpegglut-dir=${TJGDIR})
set cnfgr=($cnfgr --with-tcl-dir=${TCLDIR})
set cnfgr=($cnfgr --with-tixwish=${TIXWISH})
if ($?CPPUNITDIR) then
    set cnfgr=($cnfgr --with-cppunit-dir=${CPPUNITDIR})
endif
echo "$cnfgr" >>& $OUTPUTF
$cnfgr >>& $OUTPUTF
if ($status != 0) then
  echo "########################################################" >>& $OUTPUTF
  echo "config.log" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
  cat ${DEV_DIR}/config.log >>& $OUTPUTF
  set msg="$HOSTNAME $RELEASE_TYPE build FAILED after configure"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  rm -f ${FAILED_FILE}
  touch ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chgrp ${change_flags} fsdev ${DEV_DIR}" >>& $OUTPUTF
  chgrp ${change_flags} fsdev ${DEV_DIR} >>& $OUTPUTF
  echo "CMD: chmod ${change_flags} g+rw ${DEV_DIR}" >>& $OUTPUTF
  chmod ${change_flags} g+rw ${DEV_DIR} >>& $OUTPUTF
  chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
  chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
  exit 1
endif


#
# make
######################################################################
#
echo "##########################################################" >>& $OUTPUTF
echo "Making $DEV_DIR" >>& $OUTPUTF
echo "" >>& $OUTPUTF
echo "CMD: make -j 4" >>& $OUTPUTF
make -j 4 >>& $OUTPUTF
if ($status != 0) then
  # note: /usr/local/freesurfer/dev/bin/ dirs have not 
  # been modified (bin/ gets written after make install)
  set msg="$HOSTNAME $RELEASE_TYPE build (make) FAILED"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  rm -f ${FAILED_FILE}
  touch ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chgrp ${change_flags} fsdev ${DEV_DIR}" >>& $OUTPUTF
  chgrp ${change_flags} fsdev ${DEV_DIR} >>& $OUTPUTF
  echo "CMD: chmod ${change_flags} g+rw ${DEV_DIR}" >>& $OUTPUTF
  chmod ${change_flags} g+rw ${DEV_DIR} >>& $OUTPUTF
  chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
  chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
  exit 1  
endif


#
# make check (run available unit tests)
######################################################################
#
if ("${RELEASE_TYPE}" == "dev") then
  echo "########################################################" >>& $OUTPUTF
  echo "Make check $DEV_DIR" >>& $OUTPUTF
  echo "" >>& $OUTPUTF
#  echo "CMD: source $FREESURFER_HOME/SetUpFreeSurfer.csh" >>& $OUTPUTF
#  source $FREESURFER_HOME/SetUpFreeSurfer.csh >>& $OUTPUTF
  echo "CMD: make check" >>& $OUTPUTF
  make check >>& $OUTPUTF
  if ($status != 0) then
    # note: /usr/local/freesurfer/dev/bin/ dirs have not 
    # been modified (bin/ gets written after make install)
    set msg="$HOSTNAME $RELEASE_TYPE build (make check) FAILED unit tests"
    mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
    rm -f ${FAILED_FILE}
    touch ${FAILED_FILE}
    # set group write bit on files changed by make tools:
    echo "CMD: chgrp ${change_flags} fsdev ${DEV_DIR}" >>& $OUTPUTF
    chgrp ${change_flags} fsdev ${DEV_DIR} >>& $OUTPUTF
    echo "CMD: chmod ${change_flags} g+rw ${DEV_DIR}" >>& $OUTPUTF
    chmod ${change_flags} g+rw ${DEV_DIR} >>& $OUTPUTF
    chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
    chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
    exit 1  
  endif
endif


#
# make install
######################################################################
# (recall that configure sets $bindir to bin-new/ instead of /bin, 
# to minimize disruption of machines using contents of /bin)
echo "CMD: rm -Rf ${DEST_DIR}/bin-new" >>& $OUTPUTF
if (-e ${DEST_DIR}/bin-new) rm -rf ${DEST_DIR}/bin-new >>& $OUTPUTF
if ("${RELEASE_TYPE}" == "stable-pub") then
  # make release does make install, and runs some extra commands that
  # remove stuff not intended for public release
  echo "Building public stable" >>& $OUTPUTF
  set make_cmd=(make release)
else
  set make_cmd=(make install)
endif
echo "$make_cmd" >>& $OUTPUTF
$make_cmd >>& $OUTPUTF
if ($status != 0) then
  set msg="$HOSTNAME $RELEASE_TYPE build ($make_cmd) FAILED"
  mail -s "$msg" $FAILURE_MAIL_LIST < $OUTPUTF
  rm -f ${FAILED_FILE}
  touch ${FAILED_FILE}
  # set group write bit on files changed by make tools:
  echo "CMD: chgrp ${change_flags} fsdev ${DEV_DIR}" >>& $OUTPUTF
  chgrp ${change_flags} fsdev ${DEV_DIR} >>& $OUTPUTF
  echo "CMD: chmod ${change_flags} g+rw ${DEV_DIR}" >>& $OUTPUTF
  chmod ${change_flags} g+rw ${DEV_DIR} >>& $OUTPUTF
  chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
  chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
  # and the fsaverage in the subjects dir...
  echo "CMD: chmod ${change_flags} g+rw ${DEST_DIR}/subjects/fsaverage" >>& $OUTPUTF
  chmod ${change_flags} g+rw ${DEST_DIR}/subjects/fsaverage >>& $OUTPUTF
  chgrp ${change_flags} fsdev ${DEST_DIR}/subjects/fsaverage >>& $OUTPUTF
  exit 1  
endif
# strip symbols from binaries, greatly reducing their size
if (("${RELEASE_TYPE}" == "stable") || ("${RELEASE_TYPE}" == "stable-pub")) then
  echo "CMD: strip ${DEST_DIR}/bin-new/*" >>& $OUTPUTF
  strip ${DEST_DIR}/bin-new/* >& /dev/null
endif
#
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
echo "CMD: chgrp ${change_flags} fsdev ${DEST_DIR}" >>& $OUTPUTF
chgrp ${change_flags} fsdev ${DEST_DIR} >>& $OUTPUTF
echo "CMD: chmod ${change_flags} g+rw ${DEST_DIR}" >>& $OUTPUTF
chmod ${change_flags} g+rw ${DEST_DIR} >>& $OUTPUTF
echo "CMD: chgrp ${change_flags} fsdev ${DEV_DIR}" >>& $OUTPUTF
chgrp ${change_flags} fsdev ${DEV_DIR} >>& $OUTPUTF
echo "CMD: chmod ${change_flags} g+rw ${DEV_DIR}" >>& $OUTPUTF
chmod ${change_flags} g+rw ${DEV_DIR} >>& $OUTPUTF
chmod g+rw ${DEV_DIR}/autom4te.cache >>& $OUTPUTF
chgrp fsdev ${DEV_DIR}/config.h.in >>& $OUTPUTF
echo "CMD: chmod ${change_flags} g+rw ${LOG_DIR}" >>& $OUTPUTF
chmod ${change_flags} g+rw ${LOG_DIR} >>& $OUTPUTF


#
# library and sample subject symlinks
######################################################################
# ensure that the symlinks to the necessary packages are in place
#
symlinks:

  # first remove existing links
  rm -f ${DEST_DIR}/mni
  rm -f ${DEST_DIR}/fsl
  rm -f ${DEST_DIR}/lib/tcltktixblt
  rm -f ${DEST_DIR}/lib/gsl
  rm -f ${DEST_DIR}/lib/qt
  rm -f ${DEST_DIR}/lib/vtk
  rm -f ${DEST_DIR}/lib/vxl
  rm -f ${DEST_DIR}/lib/misc
  # then setup for proper installation
  set cmd1=(ln -s ${MNIDIR} ${DEST_DIR}/mni)
  set cmd2=(ln -s ${FSLDIR} ${DEST_DIR}/fsl)
  set cmd3=(ln -s ${TCLDIR} ${DEST_DIR}/lib/tcltktixblt)
  if ($?GSLDIR) then
    set cmd4=(ln -s ${GSLDIR} ${DEST_DIR}/lib/gsl)
  else
    set cmd4=
  endif
  if ($?VTKDIR) then
    set cmd5=(ln -s ${VTKDIR} ${DEST_DIR}/lib/vtk)
  else
    set cmd5=
  endif
  set cmd6=
  set cmd7=
  set cmd8=
  # execute the commands
  echo "$cmd1" >>& $OUTPUTF
  $cmd1
  echo "$cmd2" >>& $OUTPUTF
  $cmd2
  echo "$cmd3" >>& $OUTPUTF
  $cmd3
  echo "$cmd4" >>& $OUTPUTF
  $cmd4
  echo "$cmd5" >>& $OUTPUTF
  $cmd5
  echo "$cmd6" >>& $OUTPUTF
  $cmd6
  echo "$cmd7" >>& $OUTPUTF
  $cmd7
  echo "$cmd8" >>& $OUTPUTF
  $cmd8
  # also setup sample subject:
  rm -f ${DEST_DIR}/subjects/bert
  set cmd=(ln -s ${SPACE_FS}/subjects/bert ${DEST_DIR}/subjects/bert)
  echo "$cmd" >>& $OUTPUTF
  $cmd


#
# create build-stamp.txt
######################################################################
# create a build-stamp file, containing some basic info on this build
# which is displayed when FreeSurferEnv.csh is executed.  
# its also used by create_targz to name the tarball.
# Note: the stable build version info is hard-coded here! so it
# should be updated here with each release update
set DATE_STAMP="`date +%Y%m%d`"
set FS_PREFIX="freesurfer-${OSTYPE}-${PLATFORM}-${RELEASE_TYPE}"
if ("$RELEASE_TYPE" == "dev") then
  echo "${FS_PREFIX}-${DATE_STAMP}" > ${DEST_DIR}/build-stamp.txt
else if ("$RELEASE_TYPE" == "stable") then
  echo "${FS_PREFIX}-${STABLE_VER_NUM}-${DATE_STAMP}" \
    > ${DEST_DIR}/build-stamp.txt
else if ("$RELEASE_TYPE" == "stable-pub") then
  echo "${FS_PREFIX}-${STABLE_PUB_VER_NUM}" > ${DEST_DIR}/build-stamp.txt
else
  echo "ERROR: unknown RELEASE_TYPE: $RELEASE_TYPE"
  exit 1
endif


#
# create tarball
######################################################################
# If building stable-pub, then create a tarball
if ("$RELEASE_TYPE" == "stable-pub" || -e ${BUILD_DIR}/TARBALL ) then
  set cmd=($SCRIPT_DIR/create_targz.csh $PLATFORM $RELEASE_TYPE)
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
  if ($status) then
    echo "create_targz.csh failed to create tarball!"
    # don't exit with error, since create_targz can be re-run manually
  endif
  rm -f ${BUILD_DIR}/TARBALL >& /dev/null
endif

#
# copy libraries and include filed needed by those who want to 
# build against the freesurfer enviro (NMR center only)
######################################################################
#
if ("$RELEASE_TYPE" == "dev") then
  # remove existing 
  set cmd=(rm -Rf ${DEST_DIR}/lib/dev)
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
  set cmd=(rm -Rf ${DEST_DIR}/include)
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
  set cmd=(mkdir -p ${DEST_DIR}/lib/dev)
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
  # cp necessary libs
  set dirlist=(${DEV_DIR}/dicom/libdicom.a \
	${DEV_DIR}/hipsstubs/libhipsstubs.a \
	${DEV_DIR}/rgb/librgb.a \
	${DEV_DIR}/unix/libunix.a \
	${DEV_DIR}/utils/libutils.a)
  set cmd=(cp $dirlist ${DEST_DIR}/lib/dev)
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
  # cp include dir
  set cmd=(cp -r ${DEV_DIR}/include ${DEST_DIR})
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
  # make sure all files in DEST_DIR are group writable
  set cmd=(chmod ${change_flags} g+rw ${DEST_DIR})
  echo "$cmd" >>& $OUTPUTF
  $cmd >>& $OUTPUTF
endif


# Success, so remove fail indicator:
rm -rf ${FAILED_FILE}


done:

# Finish up
######################################################################
#
echo "##########################################################" >>& $OUTPUTF
echo "Done." >>& $OUTPUTF
set END_TIME=`date`
echo $END_TIME >>& $OUTPUTF

# Move log file to stamped version.
chmod g+w $OUTPUTF
mv $OUTPUTF ${LOG_DIR}/build_log-$RELEASE_TYPE-$HOSTNAME-$TIME_STAMP.txt
gzip -f ${LOG_DIR}/build_log-$RELEASE_TYPE-$HOSTNAME-$TIME_STAMP.txt

# Send email.
echo "Begin ${BEGIN_TIME}, end ${END_TIME}" >& $LOG_DIR/message-$HOSTNAME.txt
set msg="$HOSTNAME $RELEASE_TYPE build is wicked awesome."
mail -s "$msg" $SUCCESS_MAIL_LIST < $LOG_DIR/message-$HOSTNAME.txt
rm $LOG_DIR/message-$HOSTNAME.txt

#
# Now for a cheap way to build stable-pub, which is normally only run
# when a public distribution is needed.  Just create an empty file
# called build_stable-pub_flag in the BUILD_DIR, and stable-pub will
# be built.
if ("$RELEASE_TYPE" == "stable") then
  if (-e ${BUILD_DIR}/build_stable-pub_flag) then
    rm -f ${BUILD_DIR}/build_stable-pub_flag
    # force stable build to run again by removing a CVS'd file:
    rm -f ${DEV_DIR}/setup_configure
    ${SCRIPT_DIR}/build_stable-pub.csh
  endif
endif


# on minerva, need to clean-up some space, since /home is small
if ("$HOSTNAME" == "minerva") then
  cd ${DEV_DIR}
  make clean >>& $OUTPUTF
endif
