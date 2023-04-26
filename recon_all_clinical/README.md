# recon-all-clinical

This repository contains an implementation of recon-all-clinical, a combination of convolutional neural network and recon-all that reconstructs cortical surfaces of a clinical MRI scan of any orientation, resolution and contrast.


#### Cite the paper here and add an example figure
```
 TODO: add Karthik's preprint here!"
```
```
 B Billot, DN Greve, O Puonti, A Thielscher, K Van Leemput, B Fischl, AV Dalca, JE Iglesias:"
 SynthSeg: Domain Randomisation for Segmentation of Brain Scans of any Contrast and Resolution"
 https://arxiv.org/abs/2107.09559"
```
```
 B Billot, C Magdamo, SE Arnold, S Das, JE Iglesias:"
 Robust machine learning segmentation for large-scale analysis of heterogeneous clinical brain MRI datasets"
 https://arxiv.org/abs/2209.02032"
```
\
![Examples][def]

----------------


### General description

This tool performs recon-all-clinical, the first out-of-the-box cortical surface reconstruction and analysis of brain MRI scans of any contrast and resolution without any fine-tuning and retraining. 

This "Recon-all-like" stream for clinical scans of arbitrary orientation/resolution/contrast is essentially a combination of:
 * SynthSeg: to obtain an aseg.auto_noCCseg.mgz and to compute a Talairach transform
 * SynthSR: to have a higher resolution 1mm MPRAGE for visualization
 * SynthSurfaces: to fit surfaces by predicting the distance maps and reconstructing topologically accurate cortical surfaces

### "How to run recon-all-clinical on any clinical scan"

Once [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki) has been sourced, you can simply run **recon-all-clinical** on your own data with

```
recon-all-clinical.sh INPUT_SCAN SUBJECT_ID THREADS [SUBJECT_DIR]
```
where: 
- INPUT_SCAN: path to an image that will be processed.
- SUBJECT_ID: specifies the name or ID of the subject you would like to use. A directory with that name will be created for all the subject's FreeSurfer output.
- THREADS (optional): number of CPU threads to use. The default is just 1, so crank it up for faster processing if you have multiple cores!
- SUBJECT_DIR: only necessary if the environment variable SUBJECTS_DIR has not been set when sourcing FreeSurfer or if you want to override it.


*This stream runs a bit faster than the original recon-all, since the volumetric segmentation is much faster than the iterative Bayesian method in the standard stream*


----------------

### References

[1] *[FreeSurfer](https://www.sciencedirect.com/science/article/abs/pii/S1053811912000389?via%3Dihub)* \
Bruce Fischl \
NeuroImage (2012)



[def]: data/README_figures/examples.png