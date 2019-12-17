#!/bin/tcsh -ef

##### scripts location
setenv FSSCRIPTSDIR $FREESURFER_HOME/bin

##### IN ORDER TO ENABLE SSCNN TO RUN

##### 

##### IN ORDER TO ENABLE NIFTYREG TO RUN

setenv NIFTYREG_INSTALL /autofs/space/turan_003/users/lzollei/my_stuff/software/NiftyReg.07092019/niftyreg_install
set NiftyPath = ${NIFTYREG_INSTALL}/bin

##########################################

setenv PATH ${PATH}:${NiftyPath}:${FSSCRIPTSDIR}
echo $PATH

bash;
export PATH=/space/freesurfer/python/linux/bin:$PATH
exit
