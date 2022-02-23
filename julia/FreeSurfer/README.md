### Import this package

Before starting julia, it is recommended to define the environment variable FREESURFER_HOME. This will be used, e.g., to load the FreeSurfer color look-up table automatically.

```julia
julia> import FreeSurfer as fs
FREESURFER_HOME: /usr/local/freesurfer/dev
```

### Read .mgh, .mgz, .nii, .nii.gz volumes

```julia
julia> aa = fs.mri_read("/usr/local/freesurfer/dev/subjects/fsaverage/mri/aparc+aseg.mgz");

julia> fa = fs.mri_read("/usr/local/freesurfer/dev/trctrain/hcp/MGH35_HCP_FA_template.nii.gz");
```

### Show volume and header summary info

```julia
julia> fs.show(aa)

julia> fs.show(fa)
```

### Write .mgh, .mgz, .nii, .nii.gz volumes

```julia
julia> fs.mri_write(aa, "/tmp/aparc+aseg.nii.gz")

julia> fs.mri_write(fa, "/tmp/MGH35_HCP_FA_template.mgz")
```

### Read a .trk tractography streamline file

```julia
julia> tr = fs.trk_read("/usr/local/freesurfer/dev/trctrain/hcp/mgh_1001/syn/acomm.bbr.prep.trk");
```

### Write a .trk tractography streamline file

```julia
julia> fs.trk_write(tr, "/tmp/acomm.trk")
```

