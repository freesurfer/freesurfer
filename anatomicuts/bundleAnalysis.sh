#!/bin/sh

: ' 
Author: Alexander Zsikla
Advisor: Viviana Siless
Time: Summer 2019
Name: bundleAnalysis.sh

Description:
Takes in a patient and a possible option and produces metrics about the endpoints in order to analyze

'

# Checking correct amount of inputs
if [ $# == 0 ]; then
	echo "Usage: $0 <Patient File Directory> <Option: DTI>"
	exit 0
fi

#Making Directories
cd $1
mkdir dmri.ac
cd $1/dmri.ac
mkdir ushape
cd ushape
mkdir clusters
mkdir measures

#Saving Locations for Folders
code_location="/space/vault/7/users/vsiless/alex/Code/freesurfer/anatomicuts"
measures="$1/dmri.ac/ushape/measures"
clusters="$1/dmri.ac/ushape/clusters"

# Changing File
cd $1/mri
mri_convert wmparc.mgz wmparc.nii.gz

#Running Andrew's Code
cd $code_location
./dmri_groupByEndpoints -s $1/dmri/FOD/streamlines2anat.trk -i $1/mri/wmparc.nii.gz -d $clusters

#Running Alex's Code
if [ "$2" == "DTI" ]; then
	./dmri_extractSurfaceMeasurements -fa 4 FA $1/dmri/DTI/dti_FA.nii.gz MD $1/dmri/DTI/dti_MD.nii.gz RD $1/dmri/DTI/dti_RD.nii.gz AD $1/dmri/DTI/dti_AD.nii.gz -i $clusters/*trk -sl $1/surf/lh.white -cl $1/surf/lh.curv -tl $1/surf/lh.thickness -sr $1/surf/rh.white -cr $1/surf/rh.curv -tr $1/surf/rh.thickness -o $measures
else
	./dmri_extractSurfaceMeasurements -i $clusters/*trk -sl $1/surf/lh.white -cl $1/surf/lh.curv -tl $1/surf/lh.thickness -sr $1/surf/rh.white -cr $1/surf/rh.curv -tr $1/surf/rh.thickness -o $measures
fi

exit 0
