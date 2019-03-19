## 
##
## If You use this tool, please cite: 
## Blazejewska AI, Fischl B, Wald LL, Polimeni JR, Intracortical smoothing of small-voxel fMRI data can provide increased detection power without spatial resolution losses compared to conventional large-voxel fMRI data. NeuroImage 2019. 189:601-614.\n"
##
##

# step 0. run FreeSurfer recon-all, motion correct EPI data etc.

# step 1. generate intermediate surfaces
for l in  $(seq 0 0.1 1); do 
	lab=${l/./}
	mris_expand -thickness ${SUBJECT}/recon/surf/lh.white $l ${SUBJECT}/recon/surf/lh.midgray.$lab;
	mris_expand -thickness ${SUBJECT}/recon/surf/rh.white $l ${SUBJECT}/recon/surf/rh.midgray.$lab; 
done;

# step 2. perform registration: EPI -> STRUCTURAL


# step 3. map fMRI onto the surfaces, for each intermediate surface, for each fMRI file

MAP_SURF_DIR=${SUBJECT}/map_surf
mkdir ${MAP_SURF_DIR}

for l in ${SUBJECT}/data/ep*_mc_dt.nii.gz; do
# !!! registration matrix has to be provided
	reg=${l/%.nii.gz/_register.dat} 
	for i in $(seq 0 0.1 1); do 
		lab=${i/./}
		name=$(basename "$l"); name="${name%.*.*}"
		mri_vol2surf --mov $l --out $MAP_SURF_DIR/lh.${name}.midgray.${lab}.mgz --reg $reg --surf midgray.$lab --hemi lh
		mri_vol2surf --mov $l --out $MAP_SURF_DIR/rh.${name}.midgray.${lab}.mgz --reg $reg --surf midgray.$lab --hemi rh 
	done;
done;


# step 4. example uses of surface-based smoothing

MAP_SURF_DIR=${SUBJECT}/map_surf_smooth;
mkdir $MAP_SURF_DIR;

# with gaussian weights
mris_smooth_intracortical --surf_dir $SUBJECT/recon/surf/ --surf_name lh.midgray.05 --overlay_dir  $SUBJECT/map_surf/ --overlay_name lh.*.05.mgz --nb-rad 2 --nb-weights gauss

# with 1/NB weights (original from the manuscript)
mris_smooth_intracortical --surf_dir $SUBJECT/recon/surf/ --surf_name lh.midgray.05 --overlay_dir  $SUBJECT/map_surf/ --overlay_name lh.*.05.mgz --output_dir $MAP_SURF_DIR --output_name lh.midgray.05.nb2_1bynb.mgz --nb-rad 2 --nb-weights 1bynb

# intracortical smoothing in both direction (tangential & radial) with default gaussian smoothing
mris_smooth_intracortical --surf_dir $SUBJECT/recon/surf/ --surf_name lh.midgray.* --overlay_dir  $SUBJECT/map_surf/ --overlay_name lh.*.??.mgz --nb-rad 2 --ic-start 3 --ic-size 5 

