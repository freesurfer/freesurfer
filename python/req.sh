#!/usr/bin/env bash

func_setup_fspython()
{
   export FREESURFER_HOME=$install_path
   source $FREESURFER_HOME/SetUpFreeSurfer.sh > /dev/null 2>&1
   fspython=$install_path/bin/fspython
   if [ ! -e $fspython ]; then
      echo "$s: *** Error: Cannot find expected fspython to run as $fspython"
      exit 1
   fi
   echo -n "$s: Using fspython "
   ls $fspython
}

if [ $# == 0 ]; then
   echo "$s: Please provide one of the following arguments:"
   echo "$s: --generate   (Re)generate the requirements-build.txt and requirements-build-extra.txt files using the fspython found under your current INSTALL_PREFIX"
   echo "$s: --expand     Pass the contents of requirements file specified with -r option to the pip command line and irun the command"
   echo "$s: --add-links  Remove original requirements.txt and requirements-extra.txt and create soft links:"
   echo "$s:              requirements.txt --> requirements-build.txt"
   echo "$s:              requirements-extra.txt --> requirements-build-extra.txt"
   echo "$s: --rm-links   Remove soft links for requirements.txt and requirements-extra.txt and check out the current versions from git"
   echo "$s: --uninstall  Remove python packages that should not be re-distributed and/or are not needed"
   echo "$s: --reinstall  Re-install any python modules that were previously uninstalled (using generated postinstall.sh)"
   echo "$s: --torchcpu   Replace the existing torch package version with the cpu only version)"
   exit 0
fi

# echo "--------------------------------- start of req.sh ------------------------------"

s=`echo $0 | sed 's;^\.\/;;'`
# echo "$s: start"
generate=0
expand=0
add_links=0
rm_links=0
uninstall=0
reinstall=0
torchcpu=0

while [[ $# -gt 0 ]] && [[ "$1" == "--"* ]] ;
do
    opt="$1";
    shift;
    case "$opt" in
        "--" ) break 2;;
        "--generate" )
           echo "$s: (Re)generate the requirements-build.txt and requirements-build-extra.txt files using the fspython found under your current INSTALL_PREFIX."
           generate=1
           if [[ $expand -eq 1 || $add_links -eq 1 || $rm_links -eq 1 || $uninstall -eq 1 || $reinstall -eq 1 || $torchcpu -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        "--expand" )
           # echo "$s: Pass the contents of requirements file specified with -r option to the pip command line and run the command."
           expand=1
           if [[ $generate -eq 1 || $add_links -eq 1 || $rm_links -eq 1 || $uninstall -eq 1 || $reinstall -eq 1 || $torchcpu -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        "--add-links" )
           echo "$s: Remove original requirements.txt and requirements-extra.txt and create soft links to them from requirements-build.txt and requirements-build-extra.txt."
           add_links=1
           if [[ $expand -eq 1 || $generate -eq 1 || $rm_links -eq 1 || $uninstall -eq 1 || $reinstall -eq 1 || $torchcpu -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        "--rm-links" )
           echo "$s: Remove soft links for requirements.txt and requirements-extra.txt and check out the current versions from git."
           rm_links=1
           if [[ $expand -eq 1 || $generate -eq 1 ||$add_links -eq 1 || $uninstall -eq 1 || $reinstall -eq 1 || $torchcpu -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        "--uninstall" )
           # echo "$s: Remove packages that include libs/code we cannot re-distribute and/or are not cross-platform compatible."
           uninstall=1
           if [[ $expand -eq 1 || $generate -eq 1 || $add_links -eq 1 || $rm_links -eq 1 || $reinstall -eq 1 || $torchcpu -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        "--reinstall" )
           echo "$s: Reinstall any python modules we may have uninstalled."
           reinstall=1
           if [[ $expand -eq 1 || $generate -eq 1 || $add_links -eq 1 || $rm_links -eq 1|| $uninstall -eq 1 || $torchcpu -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        "--torchcpu" )
           echo "$s: Replace torch module with cpu only version."
           torchcpu=1
           if [[ $expand -eq 1 || $generate -eq 1 || $add_links -eq 1 || $rm_links -eq 1|| $uninstall -eq 1|| $reinstall -eq 1 ]]; then
              echo "$s: Only 1 argument allowed" && exit 1
           fi
           ;;
        *) echo >&2 "$s: Invalid option: $opt"; exit 1;;
   esac
done

this_dir=$PWD
cd ..
top_dir=$PWD
cd $this_dir

cmake_cache=""
install_path=""

if [ $expand -eq 1 ]; then
   # try to install all python packages (from various requirements files) from a single (expanded) pip command line that will be listed in the build log
   all_args="$@"
   req_file_cmd_line=`echo $all_args | sed 's;^.*-r;-r;' | awk '{print $2}'`
   remaining_cmd=`echo $all_args | sed 's;\-r[ ][ ]*'$req_file_cmd_line';;' | sed 's;[ ][ ]*$;;'`
   # Ignore any commented out package entries in the requirements file beginning with pound sign.
   #
   # Strip off +cpu from the end of any python package name if FSPYTHON_INSTALL_CUDA enabled.
   #
   # It has not turned out to be practical to make default for packages like torch = torch==2.1.2 in the python requirements file used by pip.
   # The non-cpu package version installs various gpu/NVIDIA related files (in addition to any NVIDA files installed by other packages).
   # Also the build has to work on systems w/o NVIDIA hardware/gpu's.  So the default in requirements should be, e.g., torch==2.1.2+cpui, which
   # makes less work to clean up via the --uninstall option in this script. Note that we cannot ship files with NVIDIA copyrights.

   if [ "${INTEGRATE_SAMSEG}" == "ON" ]; then
      if [ -e ${SAMSEG_REQ} ]; then
         req_files="${req_file_cmd_line} ${SAMSEG_REQ}"
      fi
   else
      req_files="${req_file_cmd_line}"
   fi

   # echo "REQ: req_files=X${req_files}X"

   if [ "${FSPYTHON_INSTALL_CUDA}" == "ON" ]; then
      echo "Stripping +cpu from the end of any package name to install since FSPYTHON_INSTALL_CUDA=${FSPYTHON_INSTALL_CUDA}"
      package_list=`cat $req_files | grep -v "^#.*" | sed 's;\+cpu$;;' | tr -s '\n' ' '`
   else
      package_list=`cat $req_files | grep -v "^#.*" | tr -s '\n' ' '`
   fi

   if [ "${FSPYTHON_INSTALL_CUDA}" == "ON" ]; then
      # does not seem to work to search this URL when listed in requirements-build-py38-cuda.txt
      str1=`echo $remaining_cmd | sed 's;pip install;pip install --extra-index-url https://download.pytorch.org/whl/cu121;'`
      expanded_cmd="$str1 $package_list"
   else 
      expanded_cmd="$remaining_cmd $package_list"
   fi
   echo "============================"
   date
   # echo "REQ FILE=${req_files}"
   # echo "REMAINING CMD=${remaining_cmd}"
   echo "ORIGINAL COMMAND:"
   echo "${all_args}"
   echo "EXPANDED WITH PACKAGES:"
   echo "${expanded_cmd}"
   echo "INSTALLING:"
   eval ${expanded_cmd}
   date
   echo "============================"
   # exit 1
fi

if [ ! -z ${BUILD_GENERATED_CMAKECACHE} ]; then
   if [ -e ${BUILD_GENERATED_CMAKECACHE} ]; then
      echo "$s: Using cmake cache file ${BUILD_GENERATED_CMAKECACHE}"
      cmake_cache=${BUILD_GENERATED_CMAKECACHE}
      grep "^CMAKE_INSTALL_PREFIX" $cmake_cache > /dev/null 2>&1
      if [ $? != 0 ]; then
         echo "$s: Cannot get install path from existing cmake cache"
         # exit 1
      else
         install_path=`grep "^CMAKE_INSTALL_PREFIX" $cmake_cache | sed 's;^.*=;;'`
      fi
   fi
fi

date_ymd=`date +%Y%m%d`
temp_install_path=""
if [ ! -z ${FS_INSTALL_DIR} ]; then
   if [ ! -e ${FS_INSTALL_DIR} ]; then
      # The nightly build install prefix has the date format above already appended, e.g., dev_20240802,
      # which will be renamed to ./dev if the build succeeds. This is a fallback in case running in
      # the nightly build tree where the build failed.
      temp_install_path="${FS_INSTALL_DIR}_${date_ymd}"
      if [ -e ${temp_install_path} ]; then
         install_path=${temp_install_path}
      fi
   else
      install_path=${FS_INSTALL_DIR}
   fi
fi

if [ "${install_path}" == "" ]; then
   echo "$s: Did not find ${FS_INSTALL_DIR}"
   echo "$s: Did not find ${temp_install_dir}"
   echo "$s: Could not determine install path for sandbox/build - exiting."
   exit 1
fi 

## Use python_binary directory as the fspython wrapper script may not yet be installed
# python_binary="${install_path}/python/bin/python3"
python_binary="$FSPYTHON_INSTALL_PREFIX/python/bin/python3"
# nvidia_subdir="$install_path/python/lib/python3.8/site-packages/nvidia"
if [ ! -z ${FSPYTHON_VERSION} ]; then
   nvidia_subdir="$FSPYTHON_INSTALL_PREFIX/python/lib/python${FSPYTHON_VERSIONi}/site-packages/nvidia"
else
   nvidia_subdir="$FSPYTHON_INSTALL_PREFIX/python/lib/python3.8/site-packages/nvidia"
fi

if [ $generate -eq 1 ]; then
   # If requirements files not soft links and modified, then stop and do not clobber
   if [ ! -L  requirements.txt ]; then
      git status -s requirements.txt | grep "M requirements.txt"
      if [ $? -eq 0 ]; then
         echo "$s: *** Error: Cannot proceed with modified requirements.txt in sandbox - please commit/push or checkout the original version."
         exit 1
      fi
   fi

   # func_setup_fspython
   build_req_new="requirements-build.txt.NEW"
   build_req_orig="requirements-build.txt.ORIG"
   build_req_git="requirements-build.txt"
   rm -f $build_req_new $build_req_orig

   # $python_binary -m pip freeze | sort | uniq > $build_req_new
   ## remove spaces around anpersand in version specs with URL's
   ## comment out entries for which pip reports no version (pyfs, qatools)
   ## comment out entries not available on MacOS (nvidia, triton)

   voxelmorph_url_when_version_invalid="voxelmorph@git+https://github.com/voxelmorph/voxelmorph.git@feb74e0541b8a390ccd2ea57b745aa8808703ca4"
   neurite_url_when_version_invalid="git+https://github.com/adalca/neurite.git@95b2b568b124cbc654467177ddcdb2cb3526788c"
   pystrum_url_when_version_invalid="git+https://github.com/adalca/pystrum.git@ba35d4b357f54e5ed577cbd413076a07ef810a21"
   spheremorph_url_when_version_invalid="spheremorph@git+https://github.com/silencer1127/spheremorph.git@master"
   ## torch cpu URL needs to be listed as arg to --find-links on line preceeding the torch cpu spec, e.g.,
   ## --find-links https://download.pytorch.org/whl/torch_stable.html
   ## torch==2.1.2+cpu
   torch_cpu_url="https://download.pytorch.org/whl/torch_stable.html"

   ## surfa now returns hash
   # surfa_url_when_version_invalid="surfa@git+https://github.com/freesurfer/surfa.git@026cabec14bb03d9dfbc6b5bdf14baec7bd51c7f"

   # $python_binary -m pip freeze | sort | uniq | sed 's; @ ;@;g' | sed 's;^qatools.*;#&;' | sed 's;^pyfs.*;#&;' | sed 's;^nvidia.*;#&;' | sed 's;^triton.*;#&;' > $build_req_new

   # $python_binary -m pip freeze | sort | uniq | sed 's; @ ;@;g' | sed 's;^qatools.*;#&;' | sed 's;^pyfs.*;#&;' | sed 's;^nvidia.*;#&;' | sed 's;^triton.*;#&;' | sed 's;voxelmorph==.*;'${voxelmorph_url_when_version_invalid}';' | sed 's;neurite==.*;'${neurite_url_when_version_invalid}';' | sed 's;pystrum==.*;'${pystrum_url_when_version_invalid}';' | sed 's;surfa==.*;'${surfa_url_when_version_invalid}';' > $build_req_new

   # $python_binary -m pip freeze | sort | uniq | sed 's; @ ;@;g' | sed 's;^qatools.*;#&;' | sed 's;^pyfs.*;#&;' | sed 's;^nvidia.*;#&;' | sed 's;^triton.*;#&;' | sed 's;voxelmorph==.*;'${voxelmorph_url_when_version_invalid}';' | sed 's;neurite==.*;'${neurite_url_when_version_invalid}';' | sed 's;pystrum==.*;'${pystrum_url_when_version_invalid}';' > $build_req_new

   ## need URL for spheremorph and addition of fsutil breaks solving requirements
   $python_binary -m pip freeze | sort | uniq | sed 's; @ ;@;g' | sed 's;^qatools.*;#&;' | sed 's;^pyfs.*;#&;' | sed 's;^nvidia.*;#&;' | sed 's;^triton.*;#&;' | sed 's;^fsutil.*;#&;' | sed 's;voxelmorph==.*;'${voxelmorph_url_when_version_invalid}';' | sed 's;neurite==.*;'${neurite_url_when_version_invalid}';' | sed 's;pystrum==.*;'${pystrum_url_when_version_invalid}';' | sed 's;spheremorph==.*;'${spheremorph_url_when_version_invalid}';' | sed 's;^torch==.*cpu;--find-links '${torch_cpu_url}'\n&;' > $build_req_new

   if [ $(wc -l < $build_req_new) -eq 0 ]; then
      echo "$s: $build_req_new has no entries so cannot use it to update requirements-build.txt"
      exit 1
   fi

   cp -p -f requirements-build.txt requirements-build.txt.ORIG
   cp -p -f $build_req_git $build_req_orig
   echo "$s: diff $build_req_orig $build_req_new"
   diff $build_req_orig $build_req_new
   # replace and let git status indicate if file modified
   cp -p -f $build_req_new $build_req_git

   cd ..

   (cd python && ls -lt $build_req_git $build_req_new $build_req_orig)
fi


if [ $add_links -eq 1 ]; then

   if [ "$cmake_cache" != "" ]; then
      grep "^FSPYTHON_BUILD_REQ:BOOL=ON" $cmake_cache > /dev/null
      if [ $? == 0 ]; then
         echo "$s: no soft links needed for requirements files with FSPYTHON_BUILD_REQ enabled"
         exit 0
      fi
   fi

   if [ -L  requirements.txt ]; then
      echo "$s: *** Error: requirements.txt is already a soft link"
      ls -l requirements.txt 
      exit 1
   fi

   if [ -L  requirements-extra.txt ]; then
      echo "$s: *** Error: requirements-extra.txt is already a soft link"
      ls -l requirements-extra.txt 
      exit 1
   fi

   if [ ! -e requirements-build.txt ]; then
      echo "$s: *** Error: requirements-build.txt does not exist to create a link from"
      exit 1
   fi

   ### NEW - zero out requirements-extra-build.txt because all revisions will be listed in requirements-build.txt
   ### file has been pushed with string below
   # rm -f requirements-extra-build.txt
   # touch requirements-extra-build.txt
   # echo "$s: # All package revisions to snapshot the python distribution are in requirements-build.txt" >> requirements-extra-build.txt

   if [ ! -e requirements-extra-build.txt ]; then
      echo "$s: *** Error: requirements-extra-build.txt does not exist to create a link from"
      exit 1
   fi

   # If requirements files modified, then stop and do not clobber
   git status -s requirements.txt | grep "M requirements.txt"
   if [ $? == 0 ]; then
      echo "$s: *** Error: Cannot proceed with modified requirements.txt in sandbox - please commit/push any changes or checkout the original version."
      exit 1
   fi
   git status -s requirements-extra.txt | grep "M requirements-extra.txt"
   if [ $? == 0 ]; then
      echo "$s: *** Error: Cannot proceed with modified requirements-extra.txt in sandbox - please commit/push any changes or checkout the original version."
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
      echo "$s: No soft link for requirements.txt to remove"
      ls -l requirements.txt 
   else
      rm -f requirements.txt
      git checkout requirements.txt
      ls -l requirements.txt
   fi
   if [ ! -L  requirements-extra.txt ]; then
      echo "$s: No soft link for requirements-extra.txt to remove"
      ls -l requirements-extra.txt 
   else
      rm -f requirements-extra.txt
      git checkout requirements-extra.txt
      ls -l requirements-extra.txt
   fi
fi

if [ $torchcpu -eq 1 ]; then
   ## Does not look like there is a way to specify the latest version of the torch cpu package
   ## in requirements-extra.txt, e.g., torch==*.cpu
   # func_setup_fspython
   torch_rev=`$python_binary -m pip freeze | grep "torch==" | sed 's;^.*==;;'`
   torch_rev_numeric=`echo $torch_rev| sed 's;\+.*;;'`
   if [ "${torch_rev}" == "${torch_rev_numeric}+cpu" ]; then
      echo "$s: ${torch_rev} already installed - nothing to do."
   else
      torch_rev_cpu="${torch_rev}+cpu"
      echo "$s: Replacing torch ${torch_rev} with torch ${torch_rev_cpu}"
      $python_binary -m pip uninstall -y torch
      if [ $? -ne 0 ]; then
         echo "$s: pip UNINSTALL failed - exiting."
         exit 1
      fi
      yes | $python_binary -m pip install torch==${torch_rev_cpu} -f https://download.pytorch.org/whl/torch_stable.html
      if [ $? -ne 0 ]; then
         echo "$s: pip INSTALL failed - exiting."
         exit 1
      fi
   fi
fi

if [ $uninstall -eq 1 ]; then
   ## check contents of nvidia directory
   # echo $s: "Contents of nvidia subdir BEFORE UNINSTALL"
   # if [ -e $nvidia_subdir ]; then ls $nvidia_subdir; fi

   rm -f postinstall.list
   # func_setup_fspython
   ## remove nvidia packages with cuda libs (installed as dependency on linux but not MacOS)
   ## remove triton
   ## replace torch with torch+cpu version (no cuda libs) via --libtorch arg above
   # $python_binary -m pip freeze | grep "^nvidia\|^triton\|^torch" > /dev/null
   $python_binary -m pip freeze | grep "^nvidia\|^triton" > /dev/null
   if [ $? -eq 0 ]; then
      if [ ! -e ./postinstall.list ]; then touch postinstall.list; fi
      # $python_binary -m pip freeze | grep '^nvidia\|^triton\|^torch' | sed 's;==.*;;' >> postinstall.list
      $python_binary -m pip freeze | grep '^nvidia\|^triton' | sed 's;==.*;;' >> postinstall.list
   else
      echo "$s: pip found nothing to uninstall for nvidia and triton in output from pip freeze."
   fi

   # tflow_subdir="$install_path/python/lib/python3.8/site-packages"
   # tflow_subdir="$install_path/python/lib/python3.8/site-packages/tensorflow"
   if [ ! -z ${FSPYTHON_VERSION} ]; then
      tflow_subdir="$FSPYTHON_INSTALL_PREFIX/python/lib/python${FSPYTHON_VERSION}/site-packages/tensorflow"
   else
      tflow_subdir="$FSPYTHON_INSTALL_PREFIX/python/lib/python3.8/site-packages/tensorflow"
   fi
   if [ -e $tflow_subdir ]; then
      rm -f header.list
      (cd ${tflow_subdir} && find -type f ! -name "*LICENSE*" ! -name "*.so*" -exec grep -i "copyright.*nvidia" {} \; -print | grep "^\.\/") > header.list
      if [ -e ./header.list ]; then
        if [ ! -z ./header.list ]; then
           file_cnt=`wc -l header.list | awk '{print $1}'`
           if [ "$file_cnt" != "0" ]; then
              echo "$s: Removing $file_cnt NVIDIA copyrighted source files"
              (cd ${tflow_subdir} && find -type f ! -name "*LICENSE*" ! -name "*.so*" -exec grep -i "copyright.*nvidia" {} \; -print | grep "^\.\/") > header.list
              rm -f delete_header.sh
              cat header.list | sed 's;^;rm -f '${tflow_subdir}'/;' > delete_header.sh
              bash delete_header.sh
              # rm -f header.list delete_header.sh
           else
              # echo "%s: no license text mentioning NVIDIA found to remove."
              echo "no license text mentioning NVIDIA found to remove."
           fi
        else
           echo "%s: no source files found to check for copyrights under $tflow_subdir"
        fi
      fi
   else
      echo "%s: tensorflow does not appear to be installed to check for NVIDIA source files"
   fi

   ## Did not work to substitute tensorflow-cpu for tensorflow, but both can be installed
   # $python_binary -m pip freeze | grep "^tensorflow==" > /dev/null
   # if [ $? -eq 0 ]; then
   #   if [ ! -e ./postinstall.list ]; then touch postinstall.list; fi
   #   $python_binary -m pip freeze | grep '^tensorflow==' | sed 's;==.*;;' >> postinstall.list
   # else
   #   echo "$s: Found nothing to uninstall for tensorflow in output from pip freeze."
   # fi

   if [ -e ./postinstall.list ]; then
      echo -n "$s: Uninstalling: "
      cat postinstall.list | tr -s '\n' ' ' && echo
      yes | $python_binary -m pip uninstall -y -q -r postinstall.list > /dev/null 2>&1
      if [ $? -ne 0 ]; then
         echo "$s: pip UNINSTALL failed - exiting."
         exit 1
      fi
   else
      # echo "%s: no postinstall.list was created."
      echo "no postinstall.list was created."
   fi

   $python_binary -m pip freeze | grep "^tensorflow" > /dev/null
   if [ $? -eq 0 ]; then
      echo "$s tensorflow modulues currently installed:"
      $python_binary -m pip freeze | grep "^tensorflow"
   else
      echo "$s: Found no tensorflow modules installed - exiting."
      exit 1
   fi

   ## check contents of nvidia directory
   # echo "$s: Contents of nvidia subdir AFTER UNINSTALL"
   # if [ -e $nvidia_subdir ]; then ls $nvidia_subdir; fi

   ## create a postinstall script to reinstall what was uninstalled
   if [ -e ./postinstall.list ]; then
      rm -f postinstall.sh
      # echo -n "$s: $python_binary -m pip install -y " >> postinstall.sh
      echo -n '$FREESURFER_HOME/python/bin/python3 -m pip install -y ' >> postinstall.sh
      # cat postinstall.list | tr -s '\n' ' ' >> postinstall.sh
      ## 03/2024 - exclude nvidia-cudnn-cu12 which breaks installation on Ubuntu linux
      cat postinstall.list | grep -v "nvidia-cudnn-cu12" | tr -s '\n' ' ' >> postinstall.sh
      chmod 755 postinstall.sh
      ## save something the fspython distribution
      cp -p -f postinstall.list $install_path/python/.
      # cp -p -f postinstall.sh $install_path/python/.
      # cat $install_path/python/postinstall.sh
      dir=`(cd . && pwd)`
      echo "$s: Created postinstall script:"
      echo "$dir/postinstall.sh"
   else
      echo "$s: postinstall.sh not created since no postinstall.list"
   fi
fi


if [ $reinstall -eq 1 ]; then
   if [ -e $install_path/python/postinstall.sh ]; then
      ## check contents of nvidia directory
      # echo "$s: Contents of nvidia subdir BEFORE REINSTALL"
      # if [ -e $nvidia_subdir ]; then ls $nvidia_subdir; fi

      func_setup_fspython
      # (cd $install_path/python && bash -x postinstall.sh)
      bash -x $install_path/python/postinstall.sh
      if [ $? -ne 0 ]; then
         echo "$s: pip REINSTALL failed - exiting."
         exit 1
      fi

      ## check contents of nvidia directory
      # echo "$s: Contents of nvidia subdir AFTER REINSTALL"
      # if [ -e $nvidia_subdir ]; then ls $nvidia_subdir; fi
   else
      echo "$s: Cannot find postinstall script postinstall.sh to re-install python modules"
   fi
fi


# echo "--------------------------------- end of req.sh ------------------------------"
# echo "$s: end"

