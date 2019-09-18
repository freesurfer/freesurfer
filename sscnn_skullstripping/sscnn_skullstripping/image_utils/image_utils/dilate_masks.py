import nibabel as nib
import glob
import os
from os.path import join as opj
import subprocess
import numpy as np
from nipype.interfaces.base import split_filename

data_dir = '/autofs/space/vault_007/users/lzollei/SkullStripping/Data'
src_subj_dir_list = sorted(glob.glob(os.path.join(data_dir, '*mprage.nii.gz')))
trg_subj_dir_list = sorted(glob.glob(os.path.join(data_dir, '*mansegmask.nii.gz')))

out_dir = '/local_mount/space/bhim/users/aj660/SkullStripping/Data'
subprocess.call(['mkdir', '-p', out_dir])

for img_file, mask_file in zip(src_subj_dir_list, trg_subj_dir_list):
    img = nib.load(img_file)
    mask = nib.load(mask_file)
    img_data = img.get_data()
    mask_data = mask.get_data()
    new_mask_data = np.zeros(mask_data.shape)


    # find extents of the msk, increase them by 32 on all sides to create a new mask
    pad_length = 32
    idx_x, idx_y, idx_z = np.where(mask_data > 0)

    min_idx_x = idx_x.min() - 32
    max_idx_x = idx_x.max() + 32



    min_idx_y = idx_y.min() - 32
    max_idx_y = idx_y.max() + 32

    min_idx_z = idx_z.min() - 32
    max_idx_z = idx_z.max() + 32


    new_mask_data[min_idx_x:max_idx_x, min_idx_y:max_idx_x, min_idx_z:max_idx_z] = 1
    img_data = img_data * new_mask_data

    out_img = nib.Nifti1Image(img_data, img.affine, img.header)
    in_file_name = split_filename(img_file)[1]
    out_file = opj(out_dir,in_file_name)
    nib.save(out_img, out_file)

    # copy the mask file as it is
    mask_file_name = split_filename(mask_file)[0]
    out_mask_file = opj(out_dir,mask_file_name)
    subprocess.call(['cp', mask_file, out_mask_file])








