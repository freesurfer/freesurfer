import nibabel as nib
import sys
import numpy as np
def filterVolumesBybval():
	if len(sys.argv)  <8:
		print("Usage {sys.argv[0]} input.nii.gz output.nii.gz bvalfileIn bvalfileOut bvecfileIn bvecfileOut  minbval maxbval")    
	print(sys.argv)

	imageIn=sys.argv[1]
	imageOut=sys.argv[2]
	bvalIn = sys.argv[3]
	bvalOut = sys.argv[4]
	bvecIn = sys.argv[5]
	bvecOut = sys.argv[6]
	minbval = float(sys.argv[7])
	maxbval = float(sys.argv[8])

	img = nib.load(imageIn)
	bs = open(bvalIn)    
	bsOut = open(bvalOut,"w") 
	bvs = open(bvecIn)    
	bvsOut = open(bvecOut,"w") 
	vecsLines= bvs.readlines()
	
	i=0
	indeces=[]   
	for row in bs:
		if float(row) >= minbval and float(row) <= maxbval:
			indeces.append(i)			
			bsOut.write(row)
			bvsOut.write(vecsLines[i])
		i+=1

	bsOut.close()
	bvsOut.close()

	print(np.shape(img.get_data()))
	data =img.get_data()[:,:,:,indeces]
	img2 = nib.Nifti1Image(data, img.get_affine())
	img2.to_filename(imageOut)

def main():
    filterVolumesBybval()
    

if __name__ == '__main__':
    main()

