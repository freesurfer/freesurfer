#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---

  DYLD_LIBRARY_PATH=".:${MCRROOT}/runtime/maci64:${MCRROOT}/bin/maci64:${MCRROOT}/sys/os/maci64:${FREESURFER_HOME}/lib/gcc/lib:${DYLD_LIBRARY_PATH}" ;
  export DYLD_LIBRARY_PATH;
  echo DYLD_LIBRARY_PATH is ${DYLD_LIBRARY_PATH};

  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} ${token}" 
      shift
  done

  RANDOMNUMBER=$(od -vAn -N4 -tu4 < /dev/urandom) ;
  MCR_CACHE_ROOT=$( echo "/tmp/MCR_${RANDOMNUMBER}/" | tr -d ' ' ) ;
  export MCR_CACHE_ROOT;
  eval "${exe_dir}/SegmentThalamicNuclei.app/Contents/MacOS/SegmentThalamicNuclei $args"
  returnVal=$?
  rm -rf $MCR_CACHE_ROOT
  
fi

exit $returnVal

