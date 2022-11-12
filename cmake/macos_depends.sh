#!/usr/bin/env bash

# STARTING WITH MACOS 12.X AND NEWER VERSIONS ...
# security now prevents linux style paths to all system libs from being visible, e.g., it is no longer
# possible to stat /usr/lib/libm.dylib. Alternate ways to reference this includes the explicit path to
# the .tbd file in the SDK distribution (see commands below) or something like the following when linking,
# -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem 
# Cmake generates correct references for system libraries like libz,
#
# mris_volmask/mris_volmask: /Library/Developer/CommandLineTools/SDKs/MacOSX12.3.sdk/usr/lib/libz.tbd
#
# - but as posts indicate, there is something "special" about a dependency on libm such that this
# does not work.  Apparently the variables that control the linker options cmake automatically adds
# to the link like "-lm" cannot be modified to change this behavior.  
#
# https://stackoverflow.com/questions/49599673/cmake-remove-added-libraries 
# https://gitlab.kitware.com/cmake/cmake/-/issues/19860

dir=$1

list_1=( $(find $dir -name "build.make") )
list_2=( $(find $dir -name "link.txt") )

for file in ${list_1[@]}
do
   grep "\/usr\/lib\/libm\.dylib" $file > /dev/null 2>&1
   if [ $? -eq 0 ]; then
       perl -i -pe's;: \/usr\/lib\/libm\.dylib;: \/Library\/Developer\/CommandLineTools\/SDKs\/MacOSX\.sdk\/usr\/lib\/libm\.tbd;g' $file
       # echo "REPLACED: /usr/lib/libm.dylib with SDK lib in $file"
   fi
done

for file in ${list_2[@]}
do
   grep "\/usr\/lib\/libm\.dylib" $file > /dev/null 2>&1
   if [ $? -eq 0 ]; then
       perl -i -pe's;\/usr\/lib\/libm\.dylib;;g' $file
       # echo "REMOVED: /usr/lib/libm.dylib from $file"
   fi
   grep " \-lm" $file > /dev/null 2>&1
   if [ $? -eq 0 ]; then
       perl -i -pe's; \-lm;;g' $file 
       # echo "REMOVED: -lm from $file"
   fi
done


