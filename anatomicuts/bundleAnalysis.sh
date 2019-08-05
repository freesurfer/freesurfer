#!/bin/sh

: ' 
Author: Alexander Zsikla
Advisor: Viviana Siless
Time: Summer 2019
Name: bundleAnalysis.sh

Description:
Takes in a directory of patients and a possible option and produces metrics about the endpoints in order to analyze

'

# Saves the location of where the code is located
code_location="/space/vault/7/users/vsiless/alex/Code/freesurfer/anatomicuts"

: '
	Name: makeDirs
	Input: The patient directory
	Returns: N/A
	Does: creates all the necessary folders within the patient directory
'
function makeDirs()
{
	cd $1
	mkdir dmri.ac
	cd dmri.ac
	mkdir ushape
	cd ushape
	mkdir clusters
	mkdir measures

	return
}

: '
	Name: run_code
	Input: a patient directory and a possible flag
	Returns: N/A
	Does: runs the code necessary for the project
'
function run_code()
{
	#Running Andrew's Code
	./dmri_groupByEndpoints -s $1/dmri/FOD/streamlines2anat.trk -i $1/mri/wmparc.nii.gz -d $clusters

	#Running Alex's Code
	if [ "$2" == "DTI" ]; then
		./dmri_extractSurfaceMeasurements -fa 4 FA $1/dmri/DTI/dti_FA2anat.nii.gz MD $1/dmri/DTI/dti_MD2anat.nii.gz RD $1/dmri/DTI/dti_RD2anat.nii.gz AD $1/dmri/DTI/dti_AD2anat.nii.gz -i $clusters/*trk -sl $1/surf/lh.pial -cl $1/surf/lh.curv -tl $1/surf/lh.thickness -sr $1/surf/rh.pial -cr $1/surf/rh.curv -tr $1/surf/rh.thickness -o $measures -ri $1/mri/wmparc.nii.gz
	else
		./dmri_extractSurfaceMeasurements -i $clusters/*trk -sl $1/surf/lh.pial -cl $1/surf/lh.curv -tl $1/surf/lh.thickness -sr $1/surf/rh.pial -cr $1/surf/rh.curv -tr $1/surf/rh.thickness -o $measures -ri $1/mri/wmparc.nii.gz
	fi

	return
}

# Checking correct amount of inputs
if [ $# == 0 ]; then
	echo "Usage: $0 <Patients Directory/*> <Option: DTI>"
	exit 0
fi

# Looks for DTI Option
DTI=
for phrase in $@; do
	if [ $phrase == "DTI" ]; then
		DTI="DTI"
	fi
done

# Cycles through all the files that are passed in
for file in $@; do
	if [ $file == "DTI" ]; then
		continue
	fi

	makeDirs $file

	# Saving Folers
	measures="$file/dmri.ac/ushape/measures"
	clusters="$file/dmri.ac/ushape/clusters"

	# Changing File
	cd $file/mri
	mri_convert wmparc.mgz wmparc.nii.gz

	# Runs the code 
	cd $code_location
	run_code $file $DTI
done

exit 0
