#!/bin/bash

#SBATCH --array=1-5
#SBATCH --time=00:28:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=4G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=BET-FAST
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023  
module load fsl/6.0.7.7 


# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_Feb.txt) 

ses_num="02"

echo "Processing subject:"
echo $sub_num

anat_folder=/home/miedemam/scratch/bids_MGH_FINAL/sub-${sub_num}/ses-${ses_num}/anat/
anat_file=${anat_folder}sub-${sub_num}_ses-${ses_num}_acq-mprage_T1w.nii.gz
out_folder=/home/miedemam/scratch/MGH_derivatives/anat/
out_file=${out_folder}sub-${sub_num}_ses-${ses_num}_acq-mprage_T1w_brain

# do the brain extraction
bet $anat_file $out_file -f 0.18 -R  

echo "Finished brain extraction!"

# also do a tissue segmentation for session 2 anatomical images only
if [ "$ses_num" = "02" ]; then
cd $out_folder
fast --channels=1 --type=1 --class=3 $out_file
echo "Finished brain segmentation!"
fi



