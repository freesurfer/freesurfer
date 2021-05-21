# SSCNN

>*Fast Infant MRI Skullstripping with Multiview 2D Convolutional Neural Networks*

### Documentation and Support
To run: 

If running on an unprocessed T1w image 

$>./sscnn_skullstrip.py  -i t1.nii.gz -o /loc/output/dir  -c t1w -gpu 0


-i: input T1w image

-o: specify output director. saved result will be /loc/output/dir/sscnn_skullstrip.nii.gz

-c: specify contrast. t1w or t2w

-gpu: for GPU to run on. Set to -1 if not running on a GPU.

For questions and support, email amoddoma@gmail.com

### Installation
Requires python3. More detailed instructions to come soon.

### License

Terms and conditions for use, reproduction, distribution and contribution are found in the FreeSurfer software license agreement contained in the **LICENSE.txt** file and [on our wiki](https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense).

*FreeSurfer © 2021 The General Hospital Corporation (Boston, MA) "MGH"*
