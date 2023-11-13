# Add this to the pip requirement and run `pip install -r requirement.txt`
spheremorph @ git+https://github.com/silencer1127/spheremorph.git@master

# Usage:
```
conda activate your_fsdev_env
cd freesurfer_josa
python fs_josa.py --hmsp lh --sub_dir /autofs/space/curv_001/users/jli/sub_test --model ./fs_josa_lh_20231025.h5 --output_dir /autofs/space/curv_001/users/jli/output_test
```

# Notes:
- This program assume `?h.sphere.rot` is avaiable in the subject/surf directory. If not, we can run the following in shell to generate
```
mris_register -n 0 ?h.sphere $FREESURFER_HOME/average/?h.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ?h.sphere.rot
```

# Issues:
- Not sure how to call fsmodule within freesurfer dev
- Current tensorflow is a cpu version in the conda env (better to have the gpu version, but will not break things)
