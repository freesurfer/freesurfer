## demo.sh: demonstration script showing how to use "mris_smooth_intracortical"
##
## If you use "mris_smooth_intracortical" in your research, please cite:
## Blazejewska AI, Fischl B, Wald LL, Polimeni JR, Intracortical smoothing of small-voxel fMRI data can provide increased detection power without spatial resolution losses compared to conventional large-voxel fMRI data. NeuroImage 2019. 189:601-614.\n"
##
##

# example subject
SUBJECT=

# step 0. run FreeSurfer recon-all, motion correct EPI data etc.
# (expect user to perform this before running demo script)


# step 1. generate intermediate surfaces

for hemi in lh rh; do
  mris_expand -a 1 -s 0.05 -n 10 -thickness ${SUBJECT}/recon/surf/${hemi}.white 1 ${SUBJECT}/recon/surf/${hemi}.midgray_   3>&1 1>&2 2>&3 
done


# step 2. perform registration: EPI -> STRUCTURAL
# (expect user to perform this before running demo script)


# step 3. map fMRI onto the surfaces, for each intermediate surface, for each fMRI file

MAP_SURF_DIR=${SUBJECT}/map_surf
mkdir ${MAP_SURF_DIR}

for l in ${SUBJECT}/data/ep*_mc_dt.nii.gz; do
# !!! registration matrix has to be provided, assume is named as follows:
	reg=${l/%.nii.gz/_register.dat} 
	for i in $(seq 0 0.1 1); do 
		lab=${i/./0}
		name=$(basename "$l"); name="${name%.*.*}"
		mri_vol2surf --mov $l --out $MAP_SURF_DIR/lh.${name}.midgray_${lab}.mgz --reg $reg --surf midgray.$lab --hemi lh
		mri_vol2surf --mov $l --out $MAP_SURF_DIR/rh.${name}.midgray_${lab}.mgz --reg $reg --surf midgray.$lab --hemi rh 
	done;
done;


# step 4. example uses of surface-based smoothing

MAP_SURF_DIR=${SUBJECT}/map_surf_smooth;
mkdir $MAP_SURF_DIR;

# with gaussian weights
mris_smooth_intracortical --surf_dir $SUBJECT/recon/surf/ --surf_name lh.midgray_005 --overlay_dir  $SUBJECT/map_surf/ --overlay_name lh.*_005.mgz --tan-size 2 --tan-weights gauss

# with 1/NB weights (as used in the 2019 article)
mris_smooth_intracortical --surf_dir $SUBJECT/recon/surf/ --surf_name lh.midgray_005 --overlay_dir  $SUBJECT/map_surf/ --overlay_name lh.*_005.mgz --output_dir $MAP_SURF_DIR --output_name lh.midgray_005.nb2_1bynb.mgz --tan-size 2 --tan-weights dstance

# intracortical smoothing in both direction (tangential & radial) with default gaussian smoothing
mris_smooth_intracortical --surf_dir $SUBJECT/recon/surf/ --surf_name lh.midgray_* --overlay_dir  $SUBJECT/map_surf/ --overlay_name lh.*_???.mgz --tan-size 2 --rad-start 3 --rad-size 5 
