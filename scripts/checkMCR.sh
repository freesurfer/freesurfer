#! /bin/tcsh -f

# Eugenio: note that when I used v80, I used to look for libdctprocess.so, 
# but not I look for libmwlaunchermain.so instead (which is not present in v80)
# Same thing for 

if ( ! -e ${FREESURFER_HOME}/MCRv84/bin/glnxa64/libmwlaunchermain.so \
  && ! -e ${FREESURFER_HOME}/MCRv84/bin/maci64/libmwlaunchermain.dylib  ) then
  echo "     "
  echo "ERROR: cannot find Matlab 2014b runtime in location:"
  echo "     "
  echo "${FREESURFER_HOME}/MCRv84"
  echo "     "
  echo "It is looking for either: "
  echo "  bin/glnxa64/libmwlaunchermain.so    (Linux 64b) or"
  echo "  bin/maci64/libmwlaunchermain.dylib (Mac 64b)"
  echo " "
  echo "The hippocampal/amygdala and brainstem modules require the (free) Matlab runtime." 
  echo "You will need to download the Matlab Compiler Runtime (MCR) for Matlab 2014b." 
  echo "To do so, please run the following command (you might need root permissions):"
  echo " "
  echo "fs_install_mcr R2014b"
  echo " "
  echo "For further details, please visit https://surfer.nmr.mgh.harvard.edu/fswiki/MatlabRuntime"
  echo " "
  exit 1
endif

exit 0



