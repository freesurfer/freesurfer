
New samseg atlas that combines the standard atlas with elements of the
charm atlas (arteries, veins, better skull/xcercsf, maybe better
sinuses) with a high myelin seg of cortex with WM crowns plus vermis
and pons belly areas plus corpus callosum. These are desgined both to
allow samseg to replace the aseg more easily and to support PETsurfer.

The new atlas as the following new segs:
34 Left-WMCrowns       (2 Cerebral WM)
66 Right-WMCrowns     (41 Cerebral WM)
126 Spinal-Cord        (16 Brain-Stem)
183 Left-Vermis-Area   (8 Cerebellum GM)
184 Right-Vermis-Area (47 Cerebellum GM) 
192 Corpus_Callosum   (2 Cerebral WM, Left but could be Right)
262 Sinus (not represented in original seg)
267 Pons-Belly-Area (16 Brain-Stem)
902 Artery         (258 Head-ExtraCerebral)
907 Other-Tissues  (258 Head-ExtraCerebral)
908 Rectus-Muscles (258 Head-ExtraCerebral)
909 Mucosa         (258 Head-ExtraCerebral)
911 Skin           (258 Head-ExtraCerebral)
914 Vein           (24 CSF)
915 Bone-Cortical   (with 916 replaces 165 Skull)
916 Bone-Cancellous (with 915 replaces 165 Skull)
930 Optic-Nerve (258 Head-ExtraCerebral)
11300 ctx_lh_high_myelin  (3 Cerebral cortex)
12300 ctx_rh_high_myelin (42 Cerebral cortex)

Missing from new:
258 Head-ExtraCerebral -- replaced with other structures
165 Skull -- replaced with other structures
80 non-WM-hypointensities -- removed from the atlas

mri_binarize \
  --replaceonly 34  2 --replaceonly 66 41 \
  --replaceonly 126 16 \
  --replaceonly 183 8 --replaceonly 184 47 \
  --replaceonly 192  2\
  --replaceonly 262 0 \
  --replaceonly 267 16 \
  --replaceonly 902 258 \
  --replaceonly 907 258 \
  --replaceonly 908 258 \
  --replaceonly 909 258 \
  --replaceonly 911 258 \
  --replaceonly 914 24 \
  --replaceonly 915 165 \
  --replaceonly 916 165 \
  --replaceonly 930 258 \
  --replaceonly 11300 3 \
  --replaceonly 12300 42 \
  --i newseg.mgz --o oldseg.mgz 


20 subjects from the 39 manually labeled buckner data set here
/autofs/cluster/fssubjects/atlases/fsv6.aseg_atlas

For each subject, samseg and charm were run, and a final segmentation
was created by merging seg_edit10 (using the cortex from the surface
analysis) with the charm seg using mkwholeheadseg.m.

The samseg registration to the samseg template was used to resample
each seg into the template space. The samseg template is more-or-less
MNI152. 
/autofs/cluster/fssubjects/atlases/fsv6.aseg_atlas/samseg-atlas/mni.22.12.31

These were passed off to Stefano C, where he copied them here
/autofs/cluster/koen/stefano/Doug_Atlas/Test9109
He created the final atlas using the same 20 subjects used in the
samseg 20Subjects_smoothing2_down2_smoothingForAffine2 in the
following way 

mesh size: 9 10 9
number of upsamples: 3
stiffness: 0.1

Schedule:
- epoch_0000_out_dir: 20 iterations with edgecollapsefactor 1.0
- epoch_0001_out_dir: 20 iterations with edgecollapsefactor 1.005
- epoch_0002_out_dir: 20 iterations with edgecollapsefactor 1.025

Selected for SAMSEG atlas directory, based on the number of nodes:
- High res: epoch_0001_out_dir CurrentMeshCollection20.gz -> 58682 nodes
- Low res: epoch_0002_out_dir CurrentMeshCollection5.gz -> 24429 nodes
- Affine:  epoch_0002_out_dir CurrentMeshCollection20.gz -> 22826 nodes

The atlas directory is in: /autofs/cluster/koen/stefano/Doug_Atlas/SAMSEG_atlas_folder

For the affine registration, Stefano grouped these superclasses:

Unknown
GlobalWM 
Veins GlobalGM Putamen Pallidum GlobalCSF
Tissues Spine Spinal-Cord

==================================================================================

Note: in the compression table, I had to swap left and right vermis
indices because I had reversed them in mkwholeheadseg.m. No need to
change the "modified" or default FS ctabs.  Changed "vermis" to
"vermis-area" in modified and default FS ctabs to indicate that this
is not a fully vetted seg

Note: changed 174 "Pons" to 267 "Pons-Belly-Area" in the compression
and "modified" ctabs to reflect that is it the belly and to indicate
that it is the "area" and not to be taken as a fully vetted
segmentation. I would have changed the name of 174, but it is used by
Eugenio's brainstem segmentation. Added 267 to default FS ctab.

==================================================================================
This atlas was run on the 19 subjects left out of the atlas and
compared to that of the
20Subjects_smoothing2_down2_smoothingForAffine2 atlas. In general, the
dice results were very close to the original atlas. 

=================================================================================
