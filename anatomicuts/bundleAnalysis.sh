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
	echo "Usage: $0 <Patients Directory> <Option: DTI>"
	exit 0
fi

code_location="/space/vault/7/users/vsiless/alex/Code/freesurfer/anatomicuts"

for file in $1; do
	#Making Directories
	cd $1/dmri.ac
	mkdir ushape
	cd ushape
	mkdir clusters
	mkdir measures

	# Saving Folers
	local measures="$1/$file/dmri.ac/ushape/measures"
	local clusters="$1/$file/dmri.ac/ushape/clusters"

	echo $measures
	echo $clusters

	# Changing File
	cd $1/$file/mri
	mri_convert wmparc.mgz wmparc.nii.gz
done

#Running Andrew's Code
cd $code_location
./dmri_groupByEndpoints -s $1/dmri/FOD/streamlines2anat.trk -i $1/mri/wmparc.nii.gz -d $clusters

#Running Alex's Code
if [ "$2" == "DTI" ]; then
	./dmri_extractSurfaceMeasurements -fa 4 FA $1/dmri/DTI/dti_FA2anat.nii.gz MD $1/dmri/DTI/dti_MD2anat.nii.gz RD $1/dmri/DTI/dti_RD2anat.nii.gz AD $1/dmri/DTI/dti_AD2anat.nii.gz -i $clusters/*trk -sl $1/surf/lh.pial -cl $1/surf/lh.curv -tl $1/surf/lh.thickness -sr $1/surf/rh.pial -cr $1/surf/rh.curv -tr $1/surf/rh.thickness -o $measures -ri $1/mri/wmparc.nii.gz
else
	./dmri_extractSurfaceMeasurements -i $clusters/*trk -sl $1/surf/lh.pial -cl $1/surf/lh.curv -tl $1/surf/lh.thickness -sr $1/surf/rh.pial -cr $1/surf/rh.curv -tr $1/surf/rh.thickness -o $measures -ri $1/mri/wmparc.nii.gz

fi

exit 0
