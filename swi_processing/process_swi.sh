#!/bin/bash
# A script to process swi images automatically from the DICOM images of Siemens/GE/Philips
#
# Usage is 
# process_swi <scanner_make> <basename> <dicom_file>
# or
# process_swi <scanner_make> <basename> <dicom_mag_file> <dicom_phase_file>
# where <scanner_make> is one of ge, siemens or philips
# where <basename> becomes the base name of the intermediate and the swi files generated
# The last 2 arguments differ because some makers put the Mag and Phase files in the same dicom file 
# and others put it in different dicom files. The user decides these.
#
# make sure nmr-dev-env is sourced

scanner_make="$1"
basename="$2"
bin_prefix=/Users/krish/clean2/dev/swi_processing

echo "---------------------------------------------------------"
echo "SWI Preprocessing..."
if [ "$1" == siemens ]; then
  ${bin_prefix}/swi_preprocess --scanner siemens --siemens_mag $3 --siemens_phase $4 --out_magnitude ${basename}_mag.nii --out_phase ${basename}_pha_wrapped.nii 
fi
if [ "$1" == ge ]; then
  ${bin_prefix}/swi_preprocess --scanner ge --ge_file $3 --out_magnitude ${basename}_mag.nii --out_phase ${basename}_pha_wrapped.nii 
fi
if [ "$1" == philips ]; then
  ${bin_prefix}/swi_preprocess --scanner philips --philips_file $3 --out_magnitude ${basename}_mag.nii --out_phase ${basename}_pha_wrapped.nii 
fi


echo "---------------------------------------------------------"
echo "Extracting brain mask and performing Phase Unwrapping..."
bet ${basename}_mag.nii ${basename}_mag1.nii -m
prelude -a ${basename}_mag.nii -p ${basename}_pha_wrapped.nii -v -n 6 -m ${basename}_mag1_mask.nii.gz -o ${basename}_pha_unwrapped.nii

echo "---------------------------------------------------------"
echo "SWI Processing..."
${bin_prefix}/swi_process --mag_file ${basename}_mag.nii --phase_file ${basename}_pha_unwrapped.nii.gz --swi_output ${basename}_6_swi.mgz --mip_level 6 

echo "---------------------------------------------------------"
echo "Display.."
tkmedit -f ${basename}_mag.nii -aux ${basename}_6_swi.mgz 
