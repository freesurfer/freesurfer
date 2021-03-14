# PSACNN

>*Pulse Sequence Adaptive Convolutional Neural Network for Fast Whole Brain Segmentation*

### Documentation and Support
To run: 
If working on unprocessed images make sure to add ROBEX and N4BiasCorrection to the PATH.

If running on an unprocessed T1w image 

$>./psacnn_segment.py  -i t1.nii.gz -o /loc/output/dir  -c t1w -gpu 0 -pre


-i: input T1w/T2w image

-o: specify output director. saved result will be /loc/output/dir/psacnnseg.nii.gz

-c: specify contrast. t1w or t2w

-gpu: for GPU to run on. Set to -1 if not running on a GPU.

-pre: Run preprocessing (transform to FS conformed space, run N4 bias correction, run ROBEX skullstripping) before running PSACNN


If working on images that are already in the conformed space, bias corrected, and skullstripped, run without the -pre flag.


For questions and support, email amoddoma@gmail.com

### Installation
Requires python3. More detailed instructions to come soon.

### License

Terms and conditions for use, reproduction, distribution and contribution are found in the FreeSurfer software license agreement contained in the **LICENSE.txt** file and [on our wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense).

*FreeSurfer © 2021 The General Hospital Corporation (Boston, MA) "MGH"*
