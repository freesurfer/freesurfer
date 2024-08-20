This is an implementation of a Bayesian segmentation method that relies on the histological atlas presented in the article:
"Next-Generation histological atlas and segmentation tool for echo "high-resolution in vivo human neuroimaging", by 
Casamitjana et al. 
-- preprint available at https://www.biorxiv.org/content/10.1101/2024.02.05.579016v1 

The code also relies on:

B Billot, DN Greve, O Puonti, A Thielscher, K Van Leemput, B Fischl, AV Dalca, and JE Iglesias. 
"SynthSeg: Segmentation of brain MRI scans of any contrast and resolution without retraining." 
Medical image analysis 86 (2023): 102789.

B Billot, C Magdamo, Y Cheng, SE Arnold, S Das, and JE Iglesias.
"Robust machine learning segmentation for large-scale analysis of heterogeneous clinical brain MRI datasets."
Proceedings of the National Academy of Sciences, 120 (2023): e2216399120.

JE Iglesias.
A ready-to-use machine learning tool for symmetric multi-modality registration of brain MRI. 
Scientific Reports, 13 (2023): 6657.

## Prerequisites:

This code requires a modern version of FreeSurfer (7.4.0 or newer), which must be sourced before running the code.

The first time you run the method, it will prompt you to download the atlas files, which are not distributed with the code.


## Usage:

To run the code, please use the script segment.sh as follows:

segment.sh INPUT_SCAN OUTPUT_DIRECTORY ATLAS_MODE GPU THREADS [BF_MODE] [GMM_MODE]

- INPUT SCAN: scan to process, in nii(.gz) or mgz format
- OUTPUT_DIRECTORY: directory where segmentations, volume files, etc will be written (more on this below).
- ATLAS_MODE: must be full (all 333 labels) or simplified (simpler brainstem protocol; recommended)
- GPU: set to 1 to use the GPU (*highly* recommended but requires a 24GB GPU!)
- THREADS: number of CPU threads to use (use -1 for all available threads)
- BF_MODE (optional): bias field mode: dct (default), polynomial, or hybrid
- GMM_MODE (optional): gaussian mixture model (GMM) model must be 1mm unless you define your own (see documentation)

Note that the first time that you run the code, you may be prompted you to download the atlas separately.

Also, Using a GPU (minimum memory: 24GB) is highly recommended. On the GPU, the code runs in about an hour (30 mins/hemisphere).
On the CPU, the running time depends heavily on the number of threads, but it can easily take over 10 hours if you do not
use many (>10) threads! Even if you use the GPU, we recommend using a bunch of CPU threads (e.g., 8) if possible, so the CPU 
parts of the algorithm run faster.

The default bias field mode (dct) uses a set of discrete cosine transform basis functions to model the bias field. The
polynomial mode uses a set of low-order 3D polynomials. The hybrid mode uses a combination of dct and polynomials.

The GMM model is crucial as it determines how different brain regions are grouped into tissue types for the purpose of 
image intensity modeling. This is specified though a set of files that should be found under data:

- data_[full/simplified]/gmm_components_[GMM_MODE].yaml: defines tissue classes and specificies the number of components of the corresponding GMM
- data_[full/simplified]/combined_aseg_labels_[GMM_MODE].yaml: defines the labels that belong to each tissue class
- data_[full/simplified]/combined_atlas_labels_[GMM_MODE].yaml: defines FreeSurfer ("aseg") labels that are used to initialize the parameters of each class.

We distribute a GMM_MODE named "1mm" that we have used in our experiments, and which is the default mode of the code. If you 
want to use your own model, you will need to create another triplet of files of your own (use the 1mm version as template).


## Output:

The output directory will contain the following files:

- bf_corrected.mgz: bias field corrected version of the input scan
- SynthSeg.mgz: SynthSeg segmentation of the scan at the whole structure level
- MNI_registration.mgz: deformation file with registration to MNI atlas (which can be found under data/mni.nii.gz)
- seg_[left/right].mgz: segmentation files (one per hemisphere).
- vols_[left/right].csv: files with volumes of the brain regions segmented by the atlas, in CSV format.
- lookup_table.txt: the lookup table to visualize seg_[left/right].mgz, for convenience
- done: this is an empty file that gets written upon successful completion of the pipeline.

You can visualize the output by CDing into the results directory and running the command:

freeview -v bf_corrected.mgz -v seg_left.mgz:colormap=lut:lut=lookup_table.txt -v seg_right.mgz:colormap=lut:lut=lookup_table.txt
 





