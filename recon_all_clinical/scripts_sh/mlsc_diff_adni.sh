#!/bin/bash
#SBATCH --account=lcnrtx
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%x_%j.out  
#SBATCH --array=1-206%50
#SBATCH --mail-user=kgopinath@mgh.harvard.edu
#SBATCH --mail-type=FAIL,TIME_LIMIT_90,END

cd /autofs/space/ballarat_001/users/kg149/recon_clinical

# Set input file paths
input_vol="/autofs/space/ballarat_001/users/kg149/recon_clinical/diffusion_ADNI2.txt"
input_subj="/autofs/space/ballarat_001/users/kg149/recon_clinical/diff_adni_subjid.txt"

# Extract variables from input files and set output directory
volume_inp=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $input_vol)
subj_inp=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $input_subj)
OUT_dir=$"/autofs/space/ballarat_002/users/kg149/diffusion_ADNI_recon/"

# Run recon-clinical script with variables
csh recon-all-clinical.sh $volume_inp $subj_inp 8 $OUT_dir