#!/usr/bin/env bash

export PATH="/space/freesurfer/develop:$PATH"

subdir=${1}

## NOT USED - FIX ME - not escpaed properly when passed in
# current_rpath=${2}
# current_rpath_fixed=`echo $current_rpath | sed 's;RIGIN;\$ORIGIN;g'`

qt_dir=${2}
qt_lib=`echo $qt_dir | sed 's;\/cmake.*;;'`

vtk_dir=${3}
vtk_lib=`echo $vtk_dir | sed 's;\/cmake.*;;'`

torch_dir=${4}
torch_lib=$torch_dir/lib

install_prefix=${5}
mkdir -p $install_prefix/bin

binary=`find . -name "freeview"`

# which patchelf 
# patchelf --print-rpath freeview
patchelf --set-rpath "${vtk_lib}:${qt_lib}:${torch_lib}" $binary
# patchelf --print-rpath freeview

# Defeat generated check for ORIGIN in RPATH via symlink
(cd $install_prefix/bin && rm -f freeview freeview_scribble_prompt)
cp -p -f $subdir/freeview $install_prefix/bin/freeview_scribble_prompt
(cd $install_prefix/bin && ln -s freeview_scribble_prompt freeview)
ls -l $install_prefix/bin/freeview

