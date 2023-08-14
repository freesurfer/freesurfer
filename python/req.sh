#!/usr/bin/env bash

if [ $# == 0 ]; then
   echo "Please provide one of the following arguments:"
   echo "--generate   (Re)generate the requirements-build.txt and requirements-build-extra.txt files using the fspython found under your current INSTALL_PREFIX"
   echo "--add-links  Remove original requirements.txt and requirements-extra.txt and create soft links:"
   echo "             requirements.txt --> requirements-build.txt"
   echo "             requirements-extra.txt --> requirements-build-extra.txt"
   echo "--rm-links   Remove soft links for requirements.txt and requirements-extra.txt and check out the current versions from git"
   exit 0
fi

echo "--------------------------------- start of req.sh ------------------------------"

add_links=0
rm_links=0
generate=0

while [[ $# -gt 0 ]] && [[ "$1" == "--"* ]] ;
do
    opt="$1";
    shift;
    case "$opt" in
        "--" ) break 2;;
        "--generate" )
           echo "(Re)generate the requirements-build.txt and requirements-build-extra.txt files using the fspython found under your current INSTALL_PREFIX:"
           generate=1
           if [[ $add_links -eq 1 || $rm_links -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        "--add-links" )
           echo "Remove original requirements.txt and requirements-extra.txt and create soft links to them from requirements-build.txt and requirements-build-extra.txt:"
           add_links=1
           if [[ $rm_links -eq 1 || $generate -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        "--rm-links" )
           echo "Remove soft links for requirements.txt and requirements-extra.txt and check out the current versions from git:"
           rm_links=1
           if [[ $add_links -eq 1 || $generate -eq 1 ]]; then
              echo "Only 1 argument allowed" && exit 1
           fi
           ;;
        *) echo >&2 "Invalid option: $opt"; exit 1;;
   esac
done


if [ $generate -eq 1 ]; then
   this_dir=$PWD
   cd ..
   top_dir=$PWD
   cd $this_dir

   # Cannot be running with FREESURFER_HOME already set
   if [ ! -z "${FREESURFER_HOME}" ]; then
      echo "*** Error: Cannot run this script with FREESURFER_HOME already set in the environment"
      exit 1
   fi
   if [ ! -e $top_dir/CMakeCache.txt ]; then
      echo "*** Error: Must run with an existing ./CMakeCache.txt file at the top of the freesurfer tree - CMakeCache.txt not found."
      exit 1
   fi

   # Use cached cmake output to get current install prefix else exit
   cmake_cache=$top_dir/CMakeCache.txt
   install_path=`grep "^CMAKE_INSTALL_PREFIX" $cmake_cache | sed 's;^.*=;;'`

   # Original requirements-extra.txt is needed for URL entries for some modules
   if [[ -L  requirements-extra.txt || -L requirements.txt ]]; then
      echo "*** Error: Please remove soft links with --rm-links and tun again with --generate"
      exit 1
      # rm -f requirements-extra.txt && git checkout requirements-extra.txt
   fi
   # If requirements files not soft links and modified, then stop and do not clobber
   if [ ! -L  requirements.txt ]; then
      git status -s requirements.txt | grep "M requirements.txt"
      if [ $? == 0 ]; then
         echo "*** Error: Cannot proceed with modified requirements.txt in sandbox - please commit/push or checkout the original version."
         exit 1
      fi
   fi
   if [ ! -L  requirements-extra.txt ]; then
      git status -s requirements-extra.txt | grep "M requirements-extra.txt"
      if [ $? == 0 ]; then
         echo "*** Error: Cannot proceed with modified requirements-extra.txt in sandbox - please commit/push or checkout the original version."
         exit 1
      fi
   fi
   if [ ! -e requirements-extra.txt ]; then
      echo "*** Error: Need an existing requirements-extra.txt file but none found."
      exit 1
   fi

   export FREESURFER_HOME=$install_path
   source $FREESURFER_HOME/SetUpFreeSurfer.sh > /dev/null 2>&1
   fspython=$install_path/bin/fspython
   if [ ! -e $fspython ]; then
      echo "*** Error: Cannot find expected fspython to run as $fspython"
      exit 1
   fi
   echo "Using fspython:"
   ls -l $fspython

   # Create and run small python script to run with fspython and list its package versions
   rm -f temp_req.py && touch temp_req.py
   echo "#! $fspython" >> temp_req.py
   echo "" >> temp_req.py
   echo "import pkg_resources" >> temp_req.py
   echo "installed_packages = [(d.project_name, d.version) for d in pkg_resources.working_set]" >> temp_req.py
   echo "print(installed_packages)" >> temp_req.py
   chmod 755 temp_req.py
   rm -f temp_req.list.all && touch temp_req.list.all

   # ./temp_req.py
   ./temp_req.py | sed 's;),;)\n;g' | sed 's;^\[;;' | sed 's;\]$;;' | sed 's;^ ;;' | sed 's;'\'';;g' | sed 's;[)(];;g' | sed 's;, ;==;' | sort -d | uniq  >> temp_req.list.all

   # filter out current entries in requirements-extra.txt
   rm -f temp_req.list.filter_extra && touch temp_req.list.filter_extra
   echo ".*h5py.*" >> temp_req.list.filter_extra
   echo ".*matplotlib.*" >> temp_req.list.filter_extra
   echo ".*opencv-python.*" >> temp_req.list.filter_extra
   echo ".*pandas.*" >> temp_req.list.filter_extra
   echo ".*scikit-image.*" >> temp_req.list.filter_extra
   echo ".*tensorflow.*" >> temp_req.list.filter_extra
   echo ".*torch.*" >> temp_req.list.filter_extra
   echo ".*transforms3d.*" >> temp_req.list.filter_extra
   echo ".*trimesh.*" >> temp_req.list.filter_extra

   rm -f temp_req.list.all.pruned && touch temp_req.list.all.pruned
   # remove extra package names from complete list
   grep -v -f temp_req.list.filter_extra temp_req.list.all > temp_req.list.all.pruned
   # new extra package list
   grep -f temp_req.list.filter_extra temp_req.list.all > temp_req.list.extra
   wc -l temp_req.list.all
   wc -l temp_req.list.all.pruned temp_req.list.extra

   update_files=1

   if [ $(wc -l < temp_req.list.all.pruned) -eq 0 ]; then
      echo "temp_req.list.all.pruned has no entries so cannot use it to update requirements-build.txt"
      update_files=0
   fi
   if [ $(wc -l < temp_req.list.extra) -eq 0 ]; then
      echo "temp_req.list.extra has no entries so cannot use it to update requirements-extra-build.txt"
      update_files=0
   fi

   if [ $update_files -eq 0 ]; then
      echo "*** Error: One or more requirements files failed to contain any entries after trying to update."
      exit 1
   fi

   rm -f requirements-build.txt requirements-extra-build.txt
   cat temp_req.list.all.pruned | sort -d | uniq > requirements-build.txt
   cat temp_req.list.extra | sort -d | uniq > requirements-extra-build.txt

   # FIX UP

   # Comment out reported version entries that do not work when pip tries to use them
   # pyfs==0.0.0
   perl -i -pe's;^pyfs==0.0.0;# pyfs==0.0.0;' requirements-build.txt
   # qatoolspython==0.9.6b0
   perl -i -pe's;^qatoolspython==0.9.6b0;# qatoolspython==0.9.6b0;' requirements-build.txt

   # Different entries for packages in requirements-extra.txt that don't have a version number
   # e.g., URLs containing git+https
   grep "git+https:" requirements-extra.txt >> requirements-extra-build.txt

   # get the latest commit hashes for these repos

   rm -rf temp_clone && mkdir temp_clone
   cd temp_clone

   # FIX ME: - Use of userid adalca and python module name is inconsistent in URLs

   git clone https://github.com/adalca/pystrum.git
   if [ ! -d pystrum ]; then echo "Error: *** Failed to checkout pystrum repo." && exit 1; fi
   pystrum_hash=`(cd pystrum && git log -1 | head -1 | awk '{print $2}')`
   echo "pystrum hash = $pystrum_hash"
   export PYSTRUM_HASH=$pystrum_hash

   git clone https://github.com/adalca/neurite.git
   if [ ! -d neurite ]; then echo "Error: *** Failed to checkout neurite repo." && exit 1; fi
   neurite_hash=`(cd neurite  && git log -1 | head -1 | awk '{print $2}')`
   echo "neurite hash = $neurite_hash"
   export NEURITE_HASH=$neurite_hash

   git clone https://github.com/voxelmorph/voxelmorph.git
   if [ ! -d voxelmorph ]; then echo "Error: *** Failed to checkout voxelmorph repo." && exit 1; fi
   voxelmorph_hash=`(cd voxelmorph  && git log -1 | head -1 | awk '{print $2}')`
   echo "voxelmorph hash = $voxelmorph_hash"
   export VOXELMORPH_HASH=$voxelmorph_hash

   cd ..

   # FIX ME: - Do a better substituttion to go from,
   #
   #  # git+https://github.com/voxelmorph/voxelmorph.git@80d0c489febfb4fa32b4a247629e79720fbb4c14
   #  voxelmorph @ git+https://github.com/voxelmorph/voxelmorph.git@dev
   #
   #  ... to ...
   #
   #  git+https://github.com/voxelmorph/voxelmorph.git@9fa9cd7631741a73a374d22d6a62ef63afdd23b4
   #  # voxelmorph @ git+https://github.com/voxelmorph/voxelmorph.git@dev

   # ucomment exsting entries starting with git+https
   perl -i -pe's;^# git\+https;git\+https;' requirements-extra-build.txt

   # comment out entries for modules that fetch latest hash and update entries with fixed hash
   perl -i -pe's;^pystrum;# pystrum;' requirements-extra-build.txt
   perl -i -pe's;^git\+https.*pystrum.*;git\+https://github.com/adalca/pystrum.git\@__PYSTRUM_HASH__;' requirements-extra-build.txt
   perl -i -pe's;__PYSTRUM_HASH__;$ENV{PYSTRUM_HASH};' requirements-extra-build.txt

   perl -i -pe's;^neurite;# neurite;' requirements-extra-build.txt
   perl -i -pe's;^git\+https.*neurite.*;git\+https://github.com/adalca/neurite.git\@__NEURITE_HASH__;' requirements-extra-build.txt
   perl -i -pe's;__NEURITE_HASH__;$ENV{NEURITE_HASH};' requirements-extra-build.txt

   perl -i -pe's;^voxelmorph;# voxelmorph;' requirements-extra-build.txt
   perl -i -pe's;^git\+https.*voxelmorph.*;git\+https://github.com/voxelmorph/voxelmorph.git\@__VOXELMORPH_HASH__;' requirements-extra-build.txt
   perl -i -pe's;__VOXELMORPH_HASH__;$ENV{VOXELMORPH_HASH};' requirements-extra-build.txt

   # diff_cmd_pystrum="diff requirements-extra.txt requirements-extra-build.txt | grep 'git+https' | grep pystrum"
   # echo && echo $diff_cmd_pystrum && eval $diff_cmd_pystrum

   # diff_cmd_neurite="diff requirements-extra.txt requirements-extra-build.txt | grep 'git+https' | grep neurite"
   # echo && echo $diff_cmd_neurite && eval $diff_cmd_neurite

   # diff_cmd_voxelmorph="diff requirements-extra.txt requirements-extra-build.txt | grep 'git+https' | grep voxelmorph"
   # echo && echo $diff_cmd_voxelmorph && eval $diff_cmd_voxelmorph
   echo

   # Generated files:
   # debug
   # ls -lt temp_req* requirements*build*
   rm -rf temp_clone
   rm -f temp_req.py temp_req.list.extra temp_req.list.all.pruned temp_req.list.filter_extra temp_req.list.all
   ls -lt requirements-build.txt requirements-extra-build.txt
fi


if [ $add_links -eq 1 ]; then
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

echo "--------------------------------- end of req.sh ------------------------------"

