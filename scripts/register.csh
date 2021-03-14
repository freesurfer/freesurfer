#! /bin/tcsh -f

#
# register.csh
#
# registers sets of COR files.
#
#
# Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
#
# Terms and conditions for use, reproduction, distribution and contribution
# are found in the 'FreeSurfer Software License Agreement' contained
# in the file 'LICENSE' found in the FreeSurfer distribution, and here:
#
# https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
#
# Reporting: freesurfer@nmr.mgh.harvard.edu
#
#


#
#
set noglob
set contents=(`cat $argv[1]/COR-.info`)
#
#

#setenv PATH ../`uname`/:$PATH

setenv DTYPE "3Db:0:0"
setenv DTYPE $DTYPE":"256":"256":"1":"$argv[1]/COR-???

rm /tmp/base+orig.HEAD  /tmp/base+orig.BRIK
to3d -prefix base  -session /tmp/ \
-xSLAB 128R-128L -ySLAB 128I-128S -zSLAB 128P-128A \
$DTYPE
echo testing
#
#
# set contents=(`cat $argv[2]/COR-.info`)
#
#
rm  /tmp/new+orig.HEAD  /tmp/new+orig.BRIK
to3d -prefix new   -session /tmp/ \
-xSLAB 128R-128L -ySLAB 128I-128S -zSLAB 128P-128A \
3Db:0:0:256:256:1:$argv[2]/COR-???
#
pushd /tmp/
#
rm base.res+orig.HEAD base.res+orig.BRIK
adwarp -dpar base+orig -apar base+orig \
 -prefix base.res -dxyz $argv[3]
#
rm new.res+orig.HEAD new.res+orig.BRIK
adwarp -dpar new+orig -apar new+orig \
 -prefix new.res -dxyz $argv[3]
#
echo done resizing
rm dset+orig.HEAD dset+orig.BRIK
3dvolreg -base 'base.res+orig[0]'  -verbose -verbose  -prefix dset \
-dfile regist.afni -maxite 32 new.res+orig

#
echo done registering
# adwarp it back...
# adwarp -dpar dset+orig -apar dset+orig \
#    -prefix reg.res -dxyz 1
set coords = `cat regist.afni`
#
#now rotate
rm reg.res+orig.HEAD reg.res+orig.BRIK
echo rotating $coords[2] $coords[3] $coords[4] $coords[5] $coords[6] $coords[7]
3drotate -prefix reg.res -rotate $coords[2]I $coords[3]R $coords[4]A \
 -ashift $coords[5]S $coords[6]L $coords[7]P new+orig

popd
#mkdir /tmp/registered
#rm /tmp/registered/*
#cp /tmp/reg.res+orig.BRIK /tmp/
#if test ! -d $argv[4]
    mkdir $argv[4]
#fi
cp /tmp/reg.res+orig.BRIK $argv[4]/
cp /tmp/reg.res+orig.HEAD $argv[4]/

mri_convert -raw 256 256 256 uchar /tmp/reg.res+orig.BRIK $argv[4]
#matlab <<EOF
#path(path,'/space/beowulf/3/users/omri/');
#AVW2COR('/tmp/reg.res+orig.BRIK','/tmp/registered');
#quit
#EOF
#cp -R /tmp/registered/ $argv[4]
#cp $argv[1]/COR-.info $argv[4]/
#then rm:
#
#rm /tmp/registered/*








