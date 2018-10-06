import nibabel as nib
import sys
def fixDsiStudioHeader():
    imageIn=sys.argv[1] #"/local_mount/space/namic/2/users/vsiless/data/hcp-mgh/mgh_1001/dtifit_FA_LPS.nii.gz"
    imageOut=sys.argv[2] #"/local_mount/space/namic/2/users/vsiless/data/hcp-mgh/mgh_1001/dtifit_FA_LPS_eye.nii.gz"

    print (imageIn, imageOut)

    img = nib.load(imageIn)
    
    img2 = nib.Nifti1Image(img.get_data(),img.get_header().get_base_affine())
    img2.update_header()
    img2.to_filename(imageOut)


def main():
    fixDsiStudioHeader()
    

if __name__ == '__main__':
    main()

