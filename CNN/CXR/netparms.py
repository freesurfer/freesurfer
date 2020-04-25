import numpy as np

wt_prefix = 'unet.cxr.mri'
wt_fname = wt_prefix + '.h5'

mri_unet_nfeatures = 30
mri_unet_depth = 4
nlabels = 2
mri_unet_feat_mult = 1.25
mri_unet_convs_per_level=8
dec_nf = [32, 32, 32, 32, 32, 16, 16]
enc_nf = [16, 32, 32, 32]

nlscale=.5
dec_nf = [int(32*nlscale), int(32*nlscale), int(32*nlscale), int(32*nlscale), int(32*nlscale), int(16*nlscale), int(16*nlscale)]
enc_nf = [int(16*nlscale), int(32*nlscale), int(32*nlscale), int(32*nlscale)]

fscale = 1
dec_nf_base = [32, 32, 32, 32, 16, 16]
enc_nf_base = [16, 32, 32]
enc_nf = [int(element * fscale) for element in enc_nf_base]
dec_nf = [int(element * fscale) for element in dec_nf_base]

enc_nf_affine = []
feature_scale = 1
for element in enc_nf:
    enc_nf_affine.append(element*feature_scale)
    feature_scale *= 2


#enc_nf_affine = [2*16, 4*32, 8*32, 8*32]
