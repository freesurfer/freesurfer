import nibabel as nib
import numpy as np
import sys
def removeHeader():
	imageIn=sys.argv[1] #"/local_mount/space/namic/2/users/vsiless/data/hcp-mgh/mgh_1001/dtifit_FA_LPS.nii.gz"
	imageOut=sys.argv[2] #"/local_mount/space/namic/2/users/vsiless/data/hcp-mgh/mgh_1001/dtifit_FA_LPS_eye.nii.gz"

	print (imageIn, imageOut)
	img = nib.load(imageIn)
	affine_orig= img.get_header().get_base_affine()

	affine = np.eye(4)
	for i in range(4):
		affine[i][i]=affine_orig[i][i];

	affine[1][1]*=-1

	affine[0][3]=7.4
	affine[1][3]=7.4
	affine[2][3]=11
	img2 = nib.Nifti1Image(img.get_data(), np.eye(4)) #affine)
	#img2.update_header()
	
	img2.to_filename(imageOut)


def main():
	removeHeader()
    

if __name__ == '__main__':
    main()

