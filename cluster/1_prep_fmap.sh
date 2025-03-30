#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=0:08:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=4G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=fmap_prep
#SBATCH --mail-type=ALL

module load StdEnv/2023  
module load fsl/6.0.7.7 


# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_22.txt) 

ses_num="02"

echo "Processing subject:"
echo $sub_num

fmap_folder=/home/miedemam/scratch/bids_MGH_FINAL/sub-${sub_num}/ses-${ses_num}/fmap/
mag_image=${fmap_folder}sub-${sub_num}_ses-${ses_num}_magnitude1.nii.gz
phase_image=${fmap_folder}sub-${sub_num}_ses-${ses_num}_phasediff.nii.gz

out_folder=/home/miedemam/scratch/MGH_derivatives/fmaps/
mag_image_bet=${out_folder}sub-${sub_num}_ses-${ses_num}_magnitude1_brain.nii.gz
out_image=${out_folder}sub-${sub_num}_ses-${ses_num}_fieldmap_rads.nii.gz

# do a brain extraction on fieldmap magnitudes
bet $mag_image $mag_image_bet -f 0.7

# copy over the non-BET image so that FEAT will be happy later
cp $mag_image ${out_folder}sub-${sub_num}_ses-${ses_num}_magnitude1.nii.gz

echo "Finished brain extraction for fieldmaps!"

# run preparation
echo_time=2.46 # in ms
fsl_prepare_fieldmap SIEMENS $phase_image $mag_image_bet $out_image $echo_time


echo "Finished preparing fieldmaps!"
