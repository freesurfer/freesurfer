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
ac_output=${ODMRI} 

function forAll()
{
 	function=$1 #i.e. tractography
	cluster_call="  $3 $4 $5 $6 $7 $8 $9" #"pbsubmit -n 1 -c " 
	
	cd ${SUBJECTS_DIR}                                                
	for s in */;       
	do                                                                                                                                                                         

		subject=${s//[\/]/}                                                                                                                                                 
		echo ${s}

		string="bash ${0} $1 ${subject} $2"                                                                                                                   
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
	echo ${subject}
	fdwi=${DMRI_DIR}/${subject}/dmri/data.nii.gz
	fbval=${DMRI_DIR}/${subject}/dmri/bvals
	fbvec=${DMRI_DIR}/${subject}/dmri/bvecs
	output=${DMRI_DIR}/${subject}/dmri/GQI/

	mkdir -p ${DMRI_DIR}/${subject}/dmri/GQI/

	diffusionUtils -f tractography -d ${fdwi} -b ${fbval} -v ${fbvec} -s ${output}

}
function getMaps()
{
	echo "getMaps"
	subject=$1
	if [[  ${2} == "DTI"  ]] || [[  ${2} == "DKI" ]] ; then

		model=$2 #DTI or DKI
		fdwi=${DMRI_DIR}/${subject}/dmri/data.nii.gz
		fbval=${DMRI_DIR}/${subject}/dmri/bvals
		fbvec=${DMRI_DIR}/${subject}/dmri/bvecs
		output=${DMRI_DIR}/${subject}/dmri/$model

		mkdir -p ${DMRI_DIR}/${subject}/dmri/$model
	
		#cd /space/erebus/2/users/vsiless/code/freesurfer/anatomicuts/
		#/space/freesurfer/python/linux/bin/python -c "import diffusionUtils;  diffusionUtils.getMaps($fdwi, $fbval, $fbvec,$output) " 
	
		diffusionUtils -f getMaps$model -d $fdwi -b $fbval -v $fbvec -s $output
	else
		echo "missing argument: dmri_ac.sh  getMaps SUBJECT_ID modelfit"
		echo "modelfit: DTI, DKI"
	fi

}
function anat2dwi()
{
	echo "anat2dwi"
	subject=$1
	mri_convert ${SUBJECTS_DIR}/${subject}/mri/brain.mgz  ${SUBJECT_DIR}/${subject}/mri/brain.nii.gz
	flirt -in ${SUBJECT_DIR}/${subject}/mri/brain.nii.gz -ref ${DMRI_DIR}/${subject}/dmri/GQI/gfa_map.nii.gz -omat ${DMRI_DIR}/${subject}/dmri/GQI/anat2dwi.mat

	#bbregister --s ${subject} --mov ${diffusion}/${subject}/gfa_map.nii.gz --reg ${diffusion}/${subject}/anat2dwi.dat  --t2 --o ${diffusion}/${subject}/brain2dwi.nii.gz
	#bbregister --s ${subject} --mov ${diffusion}/${subject}/gfa_map.nii.gz --reg ${diffusion}/${subject}/anat2dwi.dat --fslmat ${diffusion}/${subject}/anat2dwi.mat --t2 --o ${diffusion}/${subject}/brain2dwi.nii.gz
	mri_aparc2aseg --s ${subject} --labelwm --hypo-as-wm --rip-unknown --volmask --o ${SUBJECTS_DIR}/${subject}/mri/wm2009parc.mgz --ctxseg aparc.a2009s+aseg.mgz
	
	mri_vol2vol --mov  ${SUBJECTS_DIR}/${subject}/mri/wm2009parc.mgz --targ ${DMRI_DIR}/${subject}/dmri/GQI/gfa_map.nii.gz --o ${DMRI_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz --nearest --fsl ${DMRI_DIR}/${subject}/dmri/GQI/anat2dwi.mat #--reg ${diffusion}/${subject}/anat2dwi.dat
	mri_vol2vol --mov  ${SUBJECTS_DIR}/${subject}/mri/wmparc.mgz --targ ${DMRI_DIR}/${subject}/dmri/GQI/gfa_map.nii.gz --o ${DMRI_DIR}/${subject}/dmri/wmparc2dwi.nii.gz --nearest --fsl ${DMRI_DIR}/${subject}/dmri/GQI/anat2dwi.mat #--reg ${diffusion}/${subject}/anat2dwi.dat
}
function filterStreamlines()
{
	echo "filterStreamlines"
    subject=$1
    lenght=$2
	if [[ -e   ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
	then 
	    ${filtershortFibers} -i  ${DMRI_DIR}/${subject}/dmri/GQI/streamlines.trk -o  ${DMRI_DIR}/${subject}/dmri/GQI/streamlines_l${lenght}.trk -l ${lenght} -m ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz
	else
	    ${filtershortFibers} -i  ${DMRI_DIR}/${subject}/dmri/GQI/streamlines.trk -o  ${DMRI_DIR}/${subject}/dmri/GQI/streamlines_l${lenght}.trk -l ${lenght} -m ${SUBJECTS_DIR}/${subject}/dmri/wmparc2dwi.nii.gz
	fi
} 
#function PreAC()
#{
#	tractography $1
#	getMaps $1
#	anat2dwi $1
#	filterStreamlines $1
#}

function anatomiCuts()
{
    subject=$1
    lenght=$2
    mkdir -p ${ODMRI_DIR}/${subject}/dmri.ac/
    #rm -R ${output}/*
	if [[ -e   ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
	then 
	    string="${anatomiCutsBin} -s ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz -f ${SUBJECTS_DIR}/${subject}/dmri/GQI/streamlines_l${lenght}.trk -l a -c 200 -n 10 -e 500 -labels -o ${ODMRI_DIR}/${subject}/dmri.ac/"
	else
	    string="${anatomiCutsBin} -s ${SUBJECTS_DIR}/${subject}/dmri/wmparc2dwi.nii.gz -f ${SUBJECTS_DIR}/${subject}/dmri/GQI/streamlines_l${lenght}.trk -l a -c 200 -n 10 -e 500 -labels -o ${ODMRI_DIR}/${subject}/dmri.ac/"
	fi

    ${string}

}

function Hungarian()
{
	subject=$1
	targetSubject=$2	
	for c in ${clusters[@]};
	do	

		

		if [[ -e   ${SUBJECTS_DIR}/${subject}/dmri/wm2009parc2dwi.nii.gz ]] ;
		then
			si=${diffusion}/${targetSubject}/wm2009parc2dwi.nii.gz
			sj=${diffusion}/${subject}/wm2009parc2dwi.nii.gz
		else
			si=${diffusion}/${targetSubject}/wmparc2dwi.nii.gz
			sj=${diffusion}/${subject}/wmparc2dwi.nii.gz
		fi
		
		ci=${ac_output}/${targetSubject}
		cj=${ac_output}/${subject}
		mkdir -p ${cj}/match/	
		${HungarianBin} -s1 ${si} -s2 ${sj} -h1 ${ci}/ -h2 ${cj}/ -o ${cj}/match/${targetSubject}_${subject}_c${c}_hungarian.csv  -labels -hungarian -c ${c}
	done
}
function Measures()
{
	subject=$1
	anatomicuts=${ac_output}/${subject}
	for c in ${clusters[@]};
	do
		mkdir -p ${anatomicuts}/measures/

		string="${stats_ac_bin} -i ${anatomicuts}/  -n ${c} -c ${anatomicuts}/match/${targetSubject}_${subject}_c${c}_hungarian.csv -m 7 FA ${diffusion}/${subject}/dki_FA.nii   MD ${diffusion}/${subject}/dki_MD.nii   RD ${diffusion}/${subject}/dki_RD.nii   AD ${diffusion}/${subject}/dki_AD.nii   MK ${diffusion}/${subject}/dki_MK.nii   RK ${diffusion}/${subject}/dki_RK.nii   AK ${diffusion}/${subject}/dki_AK.nii   -o ${anatomicuts}/measures/${targetSubject}_${subject}_c${c}.csv"

		${string}
	done
}

function ToTarget()
{
        subject=$1
	output=${ac_output}/${subject}
        common="to${targetSubject}/"
	wmIn=${diffusion}/${subject}/wm2009parc2dwi.nii.gz
	wmOut=${diffusion}/${targetSubject}/wm2009parc2dwi.nii.gz
	
	mkdir -p ${diffusion_dir}/${subject}/dmri.ac/${common}/
	cp -R ${output}/* ${diffusion_dir}/${subject}/dmri.ac/
        
	common_clustering=${diffusion_dir}/${subject}/dmri.ac/${common}/
        images_clustering=/${diffusion_dir}/${subject}/dmri.ac/images/

	mkdir -p ${common_clustering}/trk
        mkdir -p ${common_clustering}/images
        mkdir ${images_clustering}
	mkdir -p ${diffusion_dir}/${subject}/dmri/${common}/
        #flirt -in ${subjects_dir}/${subject}/mri/brain.nii.gz -ref ${subjects_dir}/${targetSubject}/mri/brain.nii.gz -omat ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mat
	#fnirt --in=${subjects_dir}/${subject}/mri/brain.nii.gz --ref=${subjects_dir}/${targetSubject}/mri/brain.nii.gz --aff=${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mat --iout=${diffusion_dir}/${subject}/dmri/${common}/brainTo${targetSubject}.nii.gz --fout=${diffusion_dir}/${subject}/dmri/${common}/brainTo${targetSubject}_field.nii.gz
        
	fsl_reg ${diffusion}/${subject}/dki_FA.nii.gz ${diffusion}/${targetSubject}/dki_FA.nii.gz ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject} -e -FA 

	
	#mri_cvs_register --mov ${subject} --template ${targetSubject}  --noaseg # --nointensity 
	#mri_register   ${subjects_dir}/${subject}/mri/brain.nii.gz  ${subjects_dir}/${targetSubject}/mri/brain.nii.gz  ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mm3z



	cd ${diffusion_dir}/${subject}/dmri.ac/ 
        for f in *trk 
        do
                echo $f
                if [ ! -f ${images_clustering}/${f%.trk}.nii.gz ]; then
                        ${TractsToImageBin} -f ${f} -i ${wm} -e ${images_clustering}/${f%.trk}.nii.gz > /dev/null
                fi
        
                if [ ! -f  ${common_clustering}/images/${f%.trk}.nii.gz ]; then

                        dmri_trk2trk --in ${f} --out ${common_clustering}/trk/${f} --inref ${wmIn} --outref ${wmOut} --reg ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mat # --regnl ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}_warp.nii.gz #> /dev/null --reg ${diffusion_dir}/${subject}/dmri/${common}/dwiTo${targetSubject}.mat  
                        ${TractsToImageBin}  -f ${common_clustering}/trk/${f} -i ${wmOut} -e ${common_clustering}/images/${f%.trk}.nii.gz
                fi
        done
        cp HierarchicalHistory.csv ${common_clustering}/trk/


}

function average()
{       
	mkdir ${diffusion_dir}/average/dmri.ac/images
        correspondences="["
        imagesFolder="["
        outputFolder="\\\"${diffusion_dir}/average/dmri.ac/images/\\\""
        s2=${targetSubject}

        cd ${diffusion_dir}
        for v in */; 
        do
	        s=${v//[\/]/}
                echo $s
                if [  -f ${diffusion_dir}/$s/dmri.ac/HierarchicalHistory.csv ]; then
                        if [ ${#correspondences} -ge 3 ]; then 
                                correspondences=${correspondences}","
                                imagesFolder=${imagesFolder}","
                        fi
                        correspondences=${correspondences}"\\\"${diffusion_dir}/${s}/dmri.ac/match/${s2}_${s}_c200_hungarian.csv\\\""   
                        imagesFolder=${imagesFolder}"\\\"${diffusion_dir}/${s}/dmri.ac/to${targetSubject}/images/\\\""
                fi
        done    
        correspondences=${correspondences}"]"
        imagesFolder=${imagesFolder}"]"

        echo $correspondences
        echo $imagesFolder
        cd ${code}
        clusterIndeces=[$( echo `seq 0 1 199` | sed 's/ /,/g' )]        
        echo $clusterIndeces
        mkdir -p ${outputFolder}
        #correspondences, imagesFolder, outputFolder,  clusterIndeces   
        
        pbsubmit -n 1 -q max500 -c "python3 -c \"import anatomiCutsUtils;  anatomiCutsUtils.averageCorrespondingClusters($correspondences, $imagesFolder, $outputFolder,$clusterIndeces) \" "
        #python3 -c "import anatomiCutsUtils;  anatomiCutsUtils.averageCorrespondingClusters(${correspondences}, ${imagesFolder}, ${outputFolder},${clusterIndeces}) "
        
}

$@
