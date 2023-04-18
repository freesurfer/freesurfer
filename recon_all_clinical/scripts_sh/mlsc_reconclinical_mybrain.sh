#!/bin/bash
#SBATCH --account=lcnrtx
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=kgopinath@mgh.harvard.edu
#SBATCH --mail-type=FAIL,TIME_LIMIT_90,END

export FREESURFER_HOME=/usr/local/freesurfer/7.3.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh

csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CAHI21M.nii CAHI21M_T1_cln 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CARE01M.nii CARE01M_T1_cln 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CONT01M.nii CONT01M_T1_cln 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CAHI21M_FLAIR.nii CAHI21M_FLAIR 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CARE01M_FLAIR.nii CARE01M_FLAIR 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CONT01M_FLAIR.nii CONT01M_FLAIR 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CAHI21M_T2.nii CAHI21M_T2 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CARE01M_T2.nii CARE01M_T2 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/
csh recon-all-clinical.sh /autofs/space/ballarat_002/users/kg149/my_boys_mri/CONT01M_T2.nii CONT01M_T2 3 /autofs/space/ballarat_002/users/kg149/my_boys_mri/

# freeview -v CAHI21M/mri/orig.mgz -f CAHI21M/surf/lh.pial -f CAHI21M/surf/lh.white -v CAHI21M_FLAIR/mri/native.mgz -f CAHI21M_FLAIR/surf/lh.pial -f CAHI21M_FLAIR/surf/lh.white -v CAHI21M_T2/mri/native.mgz -f CAHI21M_T2/surf/lh.pial -f CAHI21M_T2/surf/lh.white
# freeview -v CARE01M/mri/orig.mgz -f CARE01M/surf/lh.pial -f CARE01M/surf/lh.white -v CARE01M_FLAIR/mri/native.mgz -f CARE01M_FLAIR/surf/lh.pial -f CARE01M_FLAIR/surf/lh.white -v CARE01M_T2/mri/native.mgz -f CARE01M_T2/surf/lh.pial -f CARE01M_T2/surf/lh.white
# freeview -v CONT01M/mri/orig.mgz -f CONT01M/surf/lh.pial -f CONT01M/surf/lh.white -v CONT01M_FLAIR/mri/native.mgz -f CONT01M_FLAIR/surf/lh.pial -f CONT01M_FLAIR/surf/lh.white -v CONT01M_T2/mri/native.mgz -f CONT01M_T2/surf/lh.pial -f CONT01M_T2/surf/lh.white