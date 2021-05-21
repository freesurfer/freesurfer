#!/usr/bin/env bash


#Required SUBJECTS_DIR with Freesurfer subjects directory.
#DMRI_DIR dmri_preproc directory if it is different from the SUBJECTS_DIR
#ODMRI_DIR AnatomiCuts output directory if different to SUBJECTS_DIR

#targetSubject target subject to find AnatomiCuts correspondences

if [[ ! ${SUBJECTS_DIR} ]]; then
	echo "ERROR: SUBJECTS_DIR not set!" 
	exit 1
fi

if [[ ! ${DMRI_DIR} ]]; then
	echo "hola"
	DMRI_DIR=${SUBJECTS_DIR}
fi

if [[ ! ${ODMRI_DIR} ]]; then
	ODMRI_DIR=${SUBJECTS_DIR}
fi

echo ${SUBJECTS_DIR}
echo ${DMRI_DIR}
echo ${ODMRI_DIR}

clusters=(200 150 100 50)
code=${FREESURFER_HOME}/bin/
filtershortFibers=${code}streamlineFilter
anatomiCutsBin=${code}dmri_AnatomiCuts
HungarianBin=${code}dmri_match 
stats_ac_bin=${code}dmri_stats_ac
TractsToImageBin=${code}trk_tools 
#forAll runAC subject DTI/DKI FOD/GQI 45 bb/- "pbsubmit_-n_1_-c_1"
function forAll()
{
 	function=$1 #i.e. tractography
	cluster_call=${6//_/ }
	echo ${cluster_call}
	
	cd ${SUBJECTS_DIR}                                                
	
	for s in */;       
	#for s in `ls -dir [a-zA-Z]*`;       
	do
		subject=${s//[\/]/}                                                                                                    
		echo ${subject}

		string="bash ${0} $1 ${subject} $2 $3 $4 $5"         
		echo ${string}
 		if [[ -z ${cluster_call} ]];
		then
			${string}
		else
			${cluster_call} "${string}"
		fi
	done      
}

function preprocessDWI()
{
	subject=$1
	folder=/space/snoke/1/public/vivros/data/nii_bval_bvec/
	output=/space/snoke/1/public/vivros/viv/$subject/
	mkdir -p $output
	mri_convert -f 0 $folder/$subject/jones_900/900.nii.gz $output/data_lowb.nii.gz
	bet $output/data_lowb.nii.gz $output/data_lowb_brain.nii.gz -m -f 0.2

	mri_concat --o $output/data.nii.gz --i $folder/$subject/jones_900/900.nii.gz $folder/$subject/jones_1800/1800.nii.gz
	indx=""
	num_frames=`mri_info --nframes $output/data.nii.gz`
	echo $num_frames
	
	for ((i=0; i< ${num_frames} ; ++i)); do indx="$indx 1"; done
	
	echo $indx > $output/index.txt
	echo "0 1 0 0.08" > $output/acqp.txt

	cat $folder/$subject/jones_900/900.bvals $folder/$subject/jones_1800/1800.bvals > $output/data.bvals
	cat $folder/$subject/jones_900/900.voxel_space.bvecs $folder/$subject/jones_1800/1800.voxel_space.bvecs > $output/data.bvecs
	

	eddy --mask=$output/data_lowb_brain_mask.nii.gz --imain=$output/data.nii.gz --bvecs=$output/data.bvecs --bvals=$output/data.bvals --out=$output/data_eddy --index=$output/index.txt --acqp=$output/acqp.txt

}

function tractography()
{
	echo "tractography"
	subject=$1
	model2=$2
	bb=$3
	echo ${subject}
	if [[ -e   ${DMRI_DIR}/${subject}/dmri/data.nii.gz ]] ;
	then
		fdwi=${DMRI_DIR}/${subject}/dmri/data.nii.gz
	else
		fdwi=${DMRI_DIR}/${subject}/dmri/dwi.nii.gz
	fi
	fbval=${DMRI_DIR}/${subject}/dmri/bvals
	fbvec=${DMRI_DIR}/${subject}/dmri/bvecs
	fxbvec=${DMRI_DIR}/${subject}/dmri/xbvecs
	
	if [[  ${2} == "GQI"  ]] ; then
	
		output=${DMRI_DIR}/${subject}/dmri/GQI/
		rm ${DMRI_DIR}/${subject}/dmri/GQI/*.trk
		mkdir -p ${DMRI_DIR}/${subject}/dmri/GQI/

		${FREESURFER_HOME}/bin/diffusionUtils -f tractography -d ${fdwi} -b ${fbval} -v ${fbvec} -s ${output}
	else

		${FREESURFER_HOME}/bin/diffusionUtils -f invert_bvecs -b ${fbval} -v ${fbvec} -a 0 -o ${fxbvec}
		mkdir  ${DMRI_DIR}/${subject}/dmri/FOD/
		mrconvert  -fslgrad ${fxbvec} ${fbval} -force ${fdwi} ${DMRI_DIR}/${subject}/dmri/FOD/dwi.mif 
		
		dwi2mask ${DMRI_DIR}/${subject}/dmri/FOD/dwi.mif ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif -force 
		dwi2response tournier ${DMRI_DIR}/${subject}/dmri/FOD/dwi.mif ${DMRI_DIR}/${subject}/dmri/FOD/response.txt -force 
		dwi2fod csd ${DMRI_DIR}/${subject}/dmri/FOD/dwi.mif ${DMRI_DIR}/${subject}/dmri/FOD/response.txt ${DMRI_DIR}/${subject}/dmri/FOD/fod.mif -mask ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif -force
		if [[  ${bb} == "-bb"  ]] ; then
			tckgen  -algorithm  SD_Stream -rk4 -downsample 10 -angle 50 -maxlen 150 -minlen 20 -cutoff .12 -seed_image ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif -mask ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif -select 300000  ${DMRI_DIR}/${subject}/dmri/FOD/fod.mif ${DMRI_DIR}/${subject}/dmri/FOD/tracts.tck -force
		else
			tckgen  -algorithm  SD_Stream -rk4 -downsample 10 -angle 40 -maxlen 250 -minlen 20 -cutoff .2 -seed_image ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif -mask ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif -select 300000  ${DMRI_DIR}/${subject}/dmri/FOD/fod.mif ${DMRI_DIR}/${subject}/dmri/FOD/tracts.tck -force
		fi
		
		#tckconvert ${DMRI_DIR}/${subject}/dmri/FOD/tracts.tck ${DMRI_DIR}/${subject}/dmri/FOD/tracks.vtk
	#-m ${SUBJECTS_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
	
		dwi2tensor  -mask ${DMRI_DIR}/${subject}/dmri/FOD/mask.mif ${DMRI_DIR}/${subject}/dmri/FOD/dwi.mif ${DMRI_DIR}/${subject}/dmri/FOD/tensor.mif -force
		tensor2metric -fa  ${DMRI_DIR}/${subject}/dmri/FOD/fa.mif ${DMRI_DIR}/${subject}/dmri/FOD/tensor.mif -force
		mrconvert ${DMRI_DIR}/${subject}/dmri/FOD/fa.mif ${DMRI_DIR}/${subject}/dmri/FOD/fa.nii.gz -force

		${FREESURFER_HOME}/bin/diffusionUtils -f TckToTrk -s ${DMRI_DIR}/${subject}/
	fi
}


function getMaps()
{
	echo "getMaps"
	subject=$1
	if [[  ${2} == "DTI"  ]] || [[  ${2} == "DKI" ]] ; then

		model=$2 #DTI or DKI
		if [[ -e   ${DMRI_DIR}/${subject}/dmri/data.nii.gz ]] ;
		then
			fdwi=${DMRI_DIR}/${subject}/dmri/data.nii.gz
		else
			fdwi=${DMRI_DIR}/${subject}/dmri/dwi.nii.gz
		fi
		fbval=${DMRI_DIR}/${subject}/dmri/bvals
		fbvec=${DMRI_DIR}/${subject}/dmri/bvecs
		output=${DMRI_DIR}/${subject}/dmri/$model

		mkdir -p ${DMRI_DIR}/${subject}/dmri/$model
	
		#cd /space/erebus/2/users/vsiless/code/freesurfer/anatomicuts/
		#/space/freesurfer/python/linux/bin/python -c "import diffusionUtils;  diffusionUtils.getMaps($fdwi, $fbval, $fbvec,$output) " 
	
		${FREESURFER_HOME}/bin/diffusionUtils -f getMaps$model -d $fdwi -b $fbval -v $fbvec -s $output
	else
		echo "missing argument: dmri_ac.sh  getMaps SUBJECT_ID modelfit"
		echo "modelfit: DTI, DKI"
	fi

}
function anat2dwi()
{
	echo "anat2dwi"
	subject=$1
	target=/FOD/fa.nii.gz
	mat=/FOD/anat2dwi.mat
	mri_convert ${SUBJECTS_DIR}/${subject}/mri/brain.mgz  ${SUBJECTS_DIR}/${subject}/mri/brain.nii.gz
	#flirt -in ${SUBJECTS_DIR}/${subject}/mri/brain.nii.gz -ref ${DMRI_DIR}/${subject}/dmri/GQI/gfa_map.nii.gz -omat ${DMRI_DIR}/${subject}/dmri/GQI/anat2dwi.mat
	flirt -in ${SUBJECTS_DIR}/${subject}/mri/brain.nii.gz -ref ${DMRI_DIR}/${subject}/dmri/${target} -omat ${DMRI_DIR}/${subject}/dmri/${mat}

	mri_aparc2aseg --s ${subject} --labelwm --hypo-as-wm --rip-unknown --volmask --o ${SUBJECTS_DIR}/${subject}/mri/wm2009parc.mgz --ctxseg aparc.a2009s+aseg.mgz
	
	#mri_vol2vol --mov  ${SUBJECTS_DIR}/${subject}/mri/wm2009parc.mgz --targ ${DMRI_DIR}/${subject}/dmri/GQI/gfa_map.nii.gz --o ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz --nearest --fsl ${DMRI_DIR}/${subject}/dmri/GQI/anat2dwi.mat #--reg ${diffusion}/${subject}/anat2dwi.dat
	mri_vol2vol --mov  ${SUBJECTS_DIR}/${subject}/mri/wm2009parc.mgz --targ ${DMRI_DIR}/${subject}/dmri/${target} --o ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz --nearest --fsl ${DMRI_DIR}/${subject}/dmri/${mat} #--reg ${diffusion}/${subject}/anat2dwi.dat
	#mri_vol2vol --mov  ${SUBJECTS_DIR}/${subject}/mri/wmparc.mgz --targ ${DMRI_DIR}/${subject}/dmri/GQI/gfa_map.nii.gz --o ${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz --nearest --fsl ${DMRI_DIR}/${subject}/dmri/GQI/anat2dwi.mat #--reg ${diffusion}/${subject}/anat2dwi.dat
	mri_vol2vol --mov  ${SUBJECTS_DIR}/${subject}/mri/wmparc.mgz --targ ${DMRI_DIR}/${subject}/dmri/${target} --o ${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz --nearest --fsl ${DMRI_DIR}/${subject}/dmri/${mat} #--reg ${diffusion}/${subject}/anat2dwi.dat
}

function dwi2anat()
{
	echo "dwi2anat"
	subject=$1
	model2=$2
	dmri_trk2trk --in ${DMRI_DIR}/${subject}/dmri/${model2}/streamlines.trk --out ${DMRI_DIR}/${subject}/dmri/${model2}/streamlines2anat.trk --reg ${SUBJECTS_DIR}/${subject}/mri/FSLREG.diff2struct.mat  --inref ${DMRI_DIR}/${subject}/dmri/${model2}/fa.nii.gz --outref ${SUBJECTS_DIR}/${subject}/mri/norm.mgz 
	
	mri_vol2vol --targ  ${SUBJECTS_DIR}/${subject}/mri/norm.mgz --mov ${DMRI_DIR}/${subject}/dmri/DTI/dti_FA.nii.gz --o ${DMRI_DIR}/${subject}/dmri/DTI/dti_FA2anat.nii.gz --fsl ${SUBJECTS_DIR}/${subject}/mri/FSLREG.diff2struct.mat 
	mri_vol2vol --targ  ${SUBJECTS_DIR}/${subject}/mri/norm.mgz --mov ${DMRI_DIR}/${subject}/dmri/DTI/dti_MD.nii.gz --o ${DMRI_DIR}/${subject}/dmri/DTI/dti_MD2anat.nii.gz --fsl ${SUBJECTS_DIR}/${subject}/mri/FSLREG.diff2struct.mat 
	mri_vol2vol --targ  ${SUBJECTS_DIR}/${subject}/mri/norm.mgz --mov ${DMRI_DIR}/${subject}/dmri/DTI/dti_RD.nii.gz --o ${DMRI_DIR}/${subject}/dmri/DTI/dti_RD2anat.nii.gz --fsl ${SUBJECTS_DIR}/${subject}/mri/FSLREG.diff2struct.mat 
	mri_vol2vol --targ  ${SUBJECTS_DIR}/${subject}/mri/norm.mgz --mov ${DMRI_DIR}/${subject}/dmri/DTI/dti_AD.nii.gz --o ${DMRI_DIR}/${subject}/dmri/DTI/dti_AD2anat.nii.gz --fsl ${SUBJECTS_DIR}/${subject}/mri/FSLREG.diff2struct.mat 


}

function filterStreamlines()
{
	echo "filterStreamlines"
    subject=$1
	model2=$2
    lenght=$3
	if [[ -e   ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
	then 
	    ${filtershortFibers} -i  ${DMRI_DIR}/${subject}/dmri/$model2/streamlines.trk -o  ${DMRI_DIR}/${subject}/dmri/$model2/streamlines_l${lenght}.trk -min ${lenght} -max 200 -m ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz
	else
	    echo ${filtershortFibers} -i  ${DMRI_DIR}/${subject}/dmri/$model2/streamlines.trk -o  ${DMRI_DIR}/${subject}/dmri/$model2/streamlines_l${lenght}.trk -min ${lenght} -max 200 -m ${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
	    ${filtershortFibers} -i  ${DMRI_DIR}/${subject}/dmri/$model2/streamlines.trk -o  ${DMRI_DIR}/${subject}/dmri/$model2/streamlines_l${lenght}.trk -min ${lenght} -max 200 -m ${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
	fi
} 
function runAC()
{
	subject=$1
	model=$2
	model2=$3
	lenght=$4
	bb=$5 

	tractography $subject $model2 $bb

	getMaps $subject DTI
	if [[ ${model} == "DKI" ]] ; then
		getMaps $1 DKI
	fi
	
	anat2dwi $subject $model2
	filterStreamlines $subject $model2 $lenght
	anatomiCuts $subject $model2 $lenght

}


function anatomiCuts()
{
    subject=$1
    model2=$2
    lenght=$3
    mkdir -p ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}
    rm  ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/*.trk
    rm  ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/*.vtk

	if [[ -e   ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
	then 
	    string="${anatomiCutsBin} -s ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz -f ${DMRI_DIR}/${subject}/dmri/$model2/streamlines_l${lenght}.trk -l a -c 200 -n 10 -e 500 -labels -o ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}"
	else
	    string="${anatomiCutsBin} -s ${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz -f ${DMRI_DIR}/${subject}/dmri/$model2/streamlines_l${lenght}.trk -l a -c 200 -n 10 -e 500 -labels -o ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}"
	fi

    ${string}

}
function preGA()
{
	subject=$1
	targetSubject=$2
	lenght=$3
	std=$4
	bb=$5	

	#Clean ${subject} ${targetSubject} ${lenght} ${std} ${bb}
	#Hungarian ${subject} ${targetSubject} ${lenght} ${std} ${bb}
	#Measures ${subject} ${targetSubject} ${lenght} ${std} ${bb}
	#ToAnat  ${subject} ${targetSubject} ${lenght} ${std} ${bb}
	SurfaceMeasures ${subject} ${targetSubject} ${lenght} ${std} ${bb}
	#ToTarget ${subject} ${targetSubject} ${lenght} ${std} ${bb}
	
}
function GA()
{
	targetSubject=$1
	lenght=$2
	std=$3
	labels_file=$4
	labels_cols=$5
	groupA=$6
	groupB=$7
	groups=$8
	thickness=$9
#	mkdir -p ${ODMRI_DIR}/average/dmri.ac/${2}/${3}/
#	rm 	 ${ODMRI_DIR}/average/dmri.ac/${2}/${3}/*
#${FREESURFER_HOME}/bin/anatomiCutsUtils -f GAL -m "DTI" -cf "${labels_file}" -cc ${labels_cols} -cta 200 -ts ${targetSubject} -s ${ODMRI_DIR} -d "," -ga $groupA -gb $groupB -l ${lenght} -std ${std} -pt ${groups} 
#${FREESURFER_HOME}/bin/anatomiCutsUtils -f GA -m "DTI" -cf "${labels_file}" -cc ${labels_cols} -cta 50:100:150:200 -ts ${targetSubject} -s ${ODMRI_DIR} -d "," -ga $groupA -gb $groupB -l ${lenght} -std ${std} -pt ${groups} 
${FREESURFER_HOME}/bin/anatomiCutsUtils -f diffAlongTheLine -m "DTI" -cf "${labels_file}" -cc ${labels_cols} -cta 50:100:150:200 -ts ${targetSubject} -s ${ODMRI_DIR} -d "," -ga $groupA -gb $groupB -l ${lenght} -std ${std} -pt ${groups} 
#${FREESURFER_HOME}/bin/anatomiCutsUtils -f thicknessPerStructure -m "DTI" -cf "${labels_file}" -cc ${labels_cols} -cta 50:100:150:200 -ts ${targetSubject} -s ${ODMRI_DIR} -d "," -ga $groupA -gb $groupB -l ${lenght} -std ${std} -pt ${groups} -t ${thickness} 
#${FREESURFER_HOME}/bin/anatomiCutsUtils -f connectivityGraph -m "DTI" -cf "${labels_file}" -cc ${labels_cols} -cta 50:100:150:200 -ts ${targetSubject} -s ${ODMRI_DIR} -d "," -ga $groupA -gb $groupB -l ${lenght} -std ${std} -pt ${groups} -t ${thickness} 

#${FREESURFER_HOME}/bin/anatomiCutsUtils -f GA -m "DKI" -cf "/space/snoke/1/public/vivros/data/demos_fullID.csv" -cc 0:6 -cta 200 -ts ${targetSubject} -s ${ODMRI_DIR} -d " " -ga 3 -gb 1 -l ${lenght} -std ${std} 
#${FREESURFER_HOME}/bin/anatomiCutsUtils -f GA -m "DKI" -cf "/space/snoke/1/public/vivros/data/demos_fullID.csv" -cc 0:6 -cta 200 -ts ${targetSubject} -s ${ODMRI_DIR} -d " " -ga 2 -gb 1 -l ${lenght} -std ${std} 
#${FREESURFER_HOME}/bin/anatomiCutsUtils -f GA -m "DKI" -cf "/space/snoke/1/public/vivros/data/demos_fullID.csv" -cc 0:6 -cta 200 -ts ${targetSubject} -s ${ODMRI_DIR} -d " " -ga 3 -gb 2 -l ${lenght} -std ${std} 

}

function Clean()
{
	subject=$1
	targetSubject=$2
	lenght=$3
	std=$4

	mkdir ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}

	${code}streamlineFilter -i ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/*trk -m ${ODMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz -d ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/ -c ${std}

        cp ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/HierarchicalHistory.csv ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/
}
function Hungarian()
{
	subject=$1
	targetSubject=$2	
	lenght=$3
	std=$4
	bb=$5

	for c in ${clusters[@]};
	do	
		if [[ -e   ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
		then
			si=${DMRI_DIR}/${targetSubject}/dmri/wm2009parc2dwi.nii.gz
			sj=${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz
		else
			si=${DMRI_DIR}/${targetSubject}/dmri/wmparc2dwi.nii.gz
			sj=${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
		fi
		
		ci=${ODMRI_DIR}/${targetSubject}/dmri.ac/${lenght}/${std}
		cj=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}
		mkdir -p ${cj}/match/	
		echo ${HungarianBin} -s1 ${si} -s2 ${sj} -h1 ${ci}/ -h2 ${cj}/ -o ${cj}/match/${targetSubject}_${subject}_c${c}_hungarian.csv  -labels -hungarian -c ${c} ${bb}
		${HungarianBin} -s1 ${si} -s2 ${sj} -h1 ${ci}/ -h2 ${cj}/ -o ${cj}/match/${targetSubject}_${subject}_c${c}_hungarian.csv  -labels -hungarian -c ${c} ${bb}
	done
}
function Measures()
{
	subject=$1
	targetSubject=$2	
	lenght=$3
	std=$4

	anatomicuts=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}
	diff=${ODMRI_DIR}/${subject}/dmri/
	for c in ${clusters[@]};
	do
		mkdir -p ${anatomicuts}/measures/

		if [[ -e   ${SUBJECTS_DIR}/${subject}/dmri/DKI/dki_RD.nii.gz ]] ;
		then 
			string="${stats_ac_bin} -i ${anatomicuts}/  -n ${c} -c ${anatomicuts}/match/${targetSubject}_${subject}_c${c}_hungarian.csv -m 7 FA ${diff}/DKI/dki_FA.nii.gz   MD ${diff}/DKI/dki_MD.nii.gz   RD ${diff}/DKI/dki_RD.nii.gz   AD ${diff}/DKI/dki_AD.nii.gz   MK ${diff}/DKI/dki_MK.nii.gz   RK ${diff}/DKI/dki_RK.nii.gz   AK ${diff}/DKI/dki_AK.nii.gz   -o ${anatomicuts}/measures/${targetSubject}_${subject}_c${c}.csv"
		else
			string="${stats_ac_bin} -i ${anatomicuts}/  -n ${c} -c ${anatomicuts}/match/${targetSubject}_${subject}_c${c}_hungarian.csv -m 4 FA ${diff}/DTI/dti_FA.nii.gz.gz   MD ${diff}/DTI/dti_MD.nii.gz   RD ${diff}/DTI/dti_RD.nii.gz   AD ${diff}/DTI/dti_AD.nii.gz  -o ${anatomicuts}/measures/${targetSubject}_${subject}_c${c}.csv"
		fi

		echo ${string}
		${string}
	done
}

function SurfaceMeasures()
{
	subject=$1
	targetSubject=$2	
	lenght=$3
	std=$4

	anatomicuts=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/toAnat/
	anatomicutsdiff=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/
	surf=${SUBJECTS_DIR}/${subject}/surf/
	annot=${SUBJECTS_DIR}/${subject}/label/
	mri=${SUBJECTS_DIR}/${subject}/mri/
	
	diff=${ODMRI_DIR}/${subject}/dmri/

 	${FREESURFER_HOME}/bin/dmri_extractSurfaceMeasurements -i ${anatomicuts}/*.trk -sl ${surf}/lh.pial -tl ${surf}/lh.thickness -cl ${surf}/lh.curv.pial -sr ${surf}/rh.pial -tr ${surf}/rh.thickness -cr ${surf}/rh.curv.pial -rid ${mri}/brain.nii.gz -ria ${mri}/brain.nii.gz  -al ${annot}/lh.aparc.annot -ar ${annot}/rh.aparc.annot -o ${anatomicutsdiff}/measures/	-p ${anatomicutsdiff}/match/${targetSubject}_${subject}_c200_hungarian.csv -fa 7 FA ${diff}/DKI/dki_FA.nii.gz   MD ${diff}/DKI/dki_MD.nii.gz   RD ${diff}/DKI/dki_RD.nii.gz   AD ${diff}/DKI/dki_AD.nii.gz   MK ${diff}/DKI/dki_MK.nii.gz   RK ${diff}/DKI/dki_RK.nii.gz   AK ${diff}/DKI/dki_AK.nii.gz   

}
function ToAnat()
{
       	subject=$1
	lenght=$3
	std=$4

        common="toAnat/"
	anatomicuts=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}
	diff=${ODMRI_DIR}/${subject}/dmri/
	surf=${SUBJECTS_DIR}/${subject}/surf/
	annot=${SUBJECTS_DIR}/${subject}/label/
	mri=${SUBJECTS_DIR}/${subject}/mri/
	

	if [[ -e   ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
	then 
		wmIn=${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz
		wmOut=${SUBJECTS_DIR}/${subject}/mri/brain.nii.gz
	else
		wmIn=${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
		wmOut=${SUBJECTS_DIR}/${subject}/mri/brain.nii.gz
	fi
	mkdir -p ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/${common}/
        
	common_clustering=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/${common}/
	mkdir ${diff}/xfms/

	flirt -in ${diff}/DTI/dti_FA.nii.gz -ref ${mri}/brain.nii.gz  -omat ${diff}/xfms/fa2brain.mat

	cd ${DMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}
        for f in *trk 
        do
                echo $f

                dmri_trk2trk --in ${f} --out ${common_clustering}/${f} --inref ${wmIn} --outref ${wmOut} --reg ${diff}/xfms/fa2brain.mat  
		/usr/pubsw/packages/dtk/current/track_info ${common_clustering}/${f} -vorder LAS LAS 
        done
        cp HierarchicalHistory.csv ${common_clustering}/

}


function ToTarget()
{
       	subject=$1
	targetSubject=$2	
	lenght=$3
	std=$4

        common="to${targetSubject}/"

	if [[ -e   ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
	then 
		wmIn=${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz
		wmOut=${DMRI_DIR}/${targetSubject}/dmri/wm2009parc2dwi.nii.gz
	else
		wmIn=${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
		wmOut=${DMRI_DIR}/${targetSubject}/dmri/wmparc2dwi.nii.gz
	fi
	mkdir -p ${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/${common}/
        
	common_clustering=${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/${common}/
        images_clustering=/${ODMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}/images/

	mkdir -p ${common_clustering}/trk
        mkdir -p ${common_clustering}/images
        mkdir ${images_clustering}
	mkdir -p ${DMRI_DIR}/${subject}/dmri/${common}/
        
	fsl_reg ${DMRI_DIR}/${subject}/dmri/DTI/dti_FA.nii.gz ${DMRI_DIR}/${targetSubject}/dmri/DTI/dti_FA.nii.gz ${DMRI_DIR}/${subject}/dmri/${common}/dwiTo${targetSubject} -e -FA 

	
	#mri_cvs_register --mov ${subject} --template ${targetSubject}  --noaseg # --nointensity 
	#mri_register   ${subjects_dir}/${subject}/mri/brain.nii.gz  ${subjects_dir}/${targetSubject}/mri/brain.nii.gz  ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mm3z

	cd ${DMRI_DIR}/${subject}/dmri.ac/${lenght}/${std}
        for f in *trk 
        do
                echo $f
                #if [ ! -f ${images_clustering}/${f%.trk}.nii.gz ]; then
                        ${TractsToImageBin} -f ${f} -i ${wmIn} -e ${images_clustering}/${f%.trk}.nii.gz #> /dev/null
                #fi
        
                #if [ ! -f  ${common_clustering}/images/${f%.trk}.nii.gz ]; then

                        dmri_trk2trk --in ${f} --out ${common_clustering}/trk/${f} --inref ${wmIn} --outref ${wmOut} --reg ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mat  --regnl ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}_warp.nii.gz --invnl #> /dev/null --reg ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mat  
                        ${TractsToImageBin}  -f ${common_clustering}/trk/${f} -i ${wmOut} -e ${common_clustering}/images/${f%.trk}.nii.gz
                #fi
        done
        cp HierarchicalHistory.csv ${common_clustering}/trk/

}

function average()
{       
	targetSubject=$1	
	lenght=$2
	std=$3
	labels_file=$4
	labels_cols=$5
	groupA=$6
	groupB=$7
	groups=$8
	thickness=$9
#	
	mkdir -p ${ODMRI_DIR}/average/dmri.ac/${lenght}/${std}/images
        #correspondences="["
        #imagesFolder="["
        outputFolder=${ODMRI_DIR}/average/dmri.ac/${lenght}/${std}/images/
        s2=${targetSubject}

        cd ${ODMRI_DIR}
        for v in */; 
        do
	        s=${v//[\/]/}
                echo $s
                if [  -f ${ODMRI_DIR}/$s/dmri.ac/${lenght}/HierarchicalHistory.csv ]; then
                        if [ ${#correspondences} -ge 3 ]; then 
                                correspondences=${correspondences},
                                imagesFolder=${imagesFolder},
                                subjectsNames=${subjectsNames},
                        fi
                        correspondences=${correspondences}${ODMRI_DIR}/${s}/dmri.ac/${lenght}/${std}/match/${s2}_${s}_c200_hungarian.csv   
                        imagesFolder=${imagesFolder}${ODMRI_DIR}/${s}/dmri.ac/${lenght}/${std}/to${targetSubject}/images/
			subjectsNames=${subjectsNames}${s}
			
                fi
        done    
        correspondences=${correspondences}
        imagesFolder=${imagesFolder}
	subjectsNames=${subjectsNames}	

        echo $correspondences
        echo $imagesFolder
        clusterIndeces=[$( echo `seq 0 1 199` | sed 's/ /,/g' )]        
        echo $clusterIndeces
        mkdir -p ${outputFolder}
        #correspondences, imagesFolder, outputFolder,  clusterIndeces   
        cd ${SUBJECTS_DIR} 
        #pbsubmit -n 1 -c "python3 -c \"import anatomiCutsUtils;  anatomiCutsUtils.averageCorrespondingClusters($correspondences, $imagesFolder, $outputFolder,$clusterIndeces) \" "
	${FREESURFER_HOME}/bin/anatomiCutsUtils -f averageCorrespondingClusters -co ${correspondences} -if ${imagesFolder} -of ${outputFolder} -in ${clusterIndeces} -sn ${subjectsNames} -cf "${labels_file}" -cc ${labels_cols} -d "," -pt ${groups} 
        #python3 -c "import anatomiCutsUtils;  anatomiCutsUtils.averageCorrespondingClusters(${correspondences}, ${imagesFolder}, ${outputFolder},${clusterIndeces}) "
        
}

$@

