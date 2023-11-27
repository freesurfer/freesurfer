#!/usr/bin/env bash


# Relative path to the build and installation directories
BUILD_DIR=../../build

# Set FREESURFER_HOME to the installation directory
export FREESURFER_HOME=$install_path

# Rest of your script...


func_setup_fspython()
{
   export FREESURFER_HOME=$install_path
   source $FREESURFER_HOME/SetUpFreeSurfer.sh > /dev/null 2>&1
   fspython=$FREESURFER_HOME/bin/fspython
   if [ ! -e $fspython ]; then
      echo "*** Error: Cannot find expected fspython to run as $fspython"
      exit 1
   fi
   echo "Using fspython:"
   ls -l $fspython
}

func_check_no_fshome()
{
   # Do not run if FREESURFER_HOME already set
   if [ ! -z "${FREESURFER_HOME}" ]; then
      echo "*** Error: Cannot run this script with FREESURFER_HOME already set in the environment"
      exit 1
   fi
   if [ ! -e $top_dir/CMakeCache.txt ]; then
      echo "*** Error: Must run with an existing ./CMakeCache.txt file at the top of the freesurfer tree - CMakeCache.txt not found."
      exit 1
   fi
}

if [ $# == 0 ]; then
   echo "Please provide one of the following arguments:"
   echo "--generate   (Re)generate the requirements-build.txt and requirements-build-extra.txt files using the fspython found under your current INSTALL_PREFIX"
   echo "--add-links  Remove original requirements.txt and requirements-extra.txt and create soft links:"
   echo "             requirements.txt --> requirements-build.txt"
   echo "             requirements-extra.txt --> requirements-build-extra.txt"
   echo "--rm-links   Remove soft links for requirements.txt and requirements-extra.txt and check out the current versions from git"
   echo "--uninstall  Remove soft python packages that should not be re-distributes and/or are not needed"
   exit 0
fi

echo "--------------------------------- start of req.sh ------------------------------"

add_links=0
rm_links=0
generate=0
uninstall=0

while [[ $# -gt 0 ]] && [[ "$1" == "--"* ]] ;
do
    opt="$1";
    shift;
    case "$opt" in
        "--" ) break 2;;
        "--generate" )
           echo "(Re)generate the requirements-build.txt and requirements-build-extra.txt files using the fspython found under your current INSTALL_PREFIX:"
           generate=1
           if [[ $add_links -eq 1 || $rm_links -eq 1 || $uninstall -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        "--add-links" )
           echo "Remove original requirements.txt and requirements-extra.txt and create soft links to them from requirements-build.txt and requirements-build-extra.txt:"
           add_links=1
           if [[ $generate -eq 1 || $rm_links -eq 1 || $uninstall -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        "--rm-links" )
           echo "Remove soft links for requirements.txt and requirements-extra.txt and check out the current versions from git:"
           rm_links=1
           if [[ $generate -eq 1 ||$add_links -eq 1 || $uninstall -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        "--uninstall" )
           echo "Remove packages that include libs/code we cannot re-distribute and/or are not cross-platform compatible:"
           uninstall=1
           if [[ $generate -eq 1 || $add_links -eq 1 || $rm_links -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        *) echo >&2 "Invalid option: $opt"; exit 1;;
   esac
done

this_dir=$PWD
cd ../../build
top_dir=$PWD
cd $this_dir

# Check if CMakeCache.txt exists in the build directory
if [ ! -e $top_dir/CMakeCache.txt ]; then
    echo "*** Error: CMakeCache.txt not found in the build directory."
    exit 1
fi

# Use cached cmake output to get current install prefix else exit
cmake_cache=$top_dir/CMakeCache.txt
install_path=`grep "^CMAKE_INSTALL_PREFIX" $cmake_cache | sed 's;^.*=;;'`

if [ $generate -eq 1 ]; then

   func_check_no_fshome
   # If requirements files not soft links and modified, then stop and do not clobber
   if [ ! -L  requirements.txt ]; then
      git status -s requirements.txt | grep "M requirements.txt"
      if [ $? -eq 0 ]; then
         echo "*** Error: Cannot proceed with modified requirements.txt in sandbox - please commit/push or checkout the original version."
         exit 1
      fi
   fi

   func_setup_fspython
   build_req_new="requirements-build.txt.NEW"
   build_req_orig="requirements-build.txt.ORIG"
   build_req_git="requirements-build.txt"
   rm -f $build_req_new $build_req_orig

   # $fspython -m pip freeze | sort | uniq > $build_req_new
   ## remove spaces around anpersand in version specs with URL's
   ## comment out entries for which pip reports no version (pyfs, qatools)
   ## comment out entries not available on MacOS (nvidia, triton)
   $fspython -m pip freeze | sort | uniq | sed 's; @ ;@;g' | sed 's;^qatools.*;#&;' | sed 's;^pyfs.*;#&;' | sed 's;^nvidia.*;#&;' | sed 's;^triton.*;#&;' > $build_req_new

   if [ $(wc -l < $build_req_new) -eq 0 ]; then
      echo "$build_req_new has no entries so cannot use it to update requirements-build.txt"
      exit 1
   fi

   cp -p -f requirements-build.txt requirements-build.txt.ORIG
   cp -p -f $build_req_git $build_req_orig
   echo "diff $build_req_orig $build_req_new"
   diff $build_req_orig $build_req_new
   # replace and let git status indicate if file modified
   cp -p -f $build_req_new $build_req_git

   cd ..

   (cd python && ls -lt $build_req_git $build_req_new $build_req_orig)
fi


if [ $add_links -eq 1 ]; then

   grep "^FSPYTHON_BUILD_REQ:BOOL=ON" $cmake_cache > /dev/null
   if [ $? == 0 ]; then
      echo "no soft links needed for requirements files with FSPYTHON_BUILD_REQ enabled"
      exit 0
   fi

   if [ -L  requirements.txt ]; then
      echo "*** Error: requirements.txt is already a soft link"
      ls -l requirements.txt 
      exit 1
   fi

   if [ -L  requirements-extra.txt ]; then
      echo "*** Error: requirements-extra.txt is already a soft link"
      ls -l requirements-extra.txt 
      exit 1
   fi

   if [ ! -e requirements-build.txt ]; then
      echo "*** Error: requirements-build.txt does not exist to create a link from"
      exit 1
   fi

   ### NEW - zero out requirements-extra-build.txt because all revisions will be listed in requirements-build.txt
   ### file has been pushed with string below
   # rm -f requirements-extra-build.txt
   # touch requirements-extra-build.txt
   # echo "# All package revisions to snapshot the python distribution are in requirements-build.txt" >> requirements-extra-build.txt

   if [ ! -e requirements-extra-build.txt ]; then
      echo "*** Error: requirements-extra-build.txt does not exist to create a link from"
      exit 1
   fi

   # If requirements files modified, then stop and do not clobber
   git status -s requirements.txt | grep "M requirements.txt"
   if [ $? == 0 ]; then
      echo "*** Error: Cannot proceed with modified requirements.txt in sandbox - please commit/push any changes or checkout the original version."
      exit 1
   fi
   git status -s requirements-extra.txt | grep "M requirements-extra.txt"
   if [ $? == 0 ]; then
      echo "*** Error: Cannot proceed with modified requirements-extra.txt in sandbox - please commit/push any changes or checkout the original version."
      exit 1
   fi
   rm -f requirements.txt
   ln -s requirements-build.txt requirements.txt
   ls -l requirements.txt
   rm -f requirements-extra.txt
   ln -s requirements-extra-build.txt requirements-extra.txt
   ls -l requirements-extra.txt
fi


if [ $rm_links -eq 1 ]; then
   if [ ! -L  requirements.txt ]; then
      echo "No soft link for requirements.txt to remove"
      ls -l requirements.txt 
   else
      rm -f requirements.txt
      git checkout requirements.txt
      ls -l requirements.txt
   fi
   if [ ! -L  requirements-extra.txt ]; then
      echo "No soft link for requirements-extra.txt to remove"
      ls -l requirements-extra.txt 
   else
      rm -f requirements-extra.txt
      git checkout requirements-extra.txt
      ls -l requirements-extra.txt
   fi
fi

if [ $uninstall -eq 1 ]; then

   # func_check_no_fshome
   func_setup_fspython

   # remove nvidia packages with compiled cuda shared libs (installed as dependency on linux but not MacOS)
   # remove triton
   fspython -m pip freeze | grep "^nvidia" > /dev/null
   if [ $? -eq 0 ]; then
      rm -f uninstall.txt
      fspython -m pip freeze | grep '^nvidia\|^triton' | sed 's;==.*;;' > uninstall.txt
      yes | fspython -m pip uninstall -r uninstall.txt
   else
      echo "Found nothing to uninstall in output from pip freeze."
   fi
fi

echo "--------------------------------- end of req.sh ------------------------------"

