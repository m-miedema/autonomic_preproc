#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=0:45:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=16G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=clean-func-MYSES-MYTASK
#SBATCH --mail-type=ALL
#SBATCH --output=out_clean_MYSES-MYTASK_%j.txt
module load StdEnv/2023  
module load fsl/6.0.7.7 

# choose session and task
ses_num=MYSES
task=MYTASK

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_22.txt) 

echo "Processing subject:"
in_root=sub-${sub_num}_ses-${ses_num}_${task}
bids_root=sub-${sub_num}_ses-${ses_num}_task-${task}
echo $bids_root

# do minimal cleaning
clean=minimal
Text2Vest ./regressors/${bids_root}_regressors_${clean}.txt ./regressors/fsl_des/${bids_root}_regressors_${clean}_des
# generate the cleaned data
fsl_glm -i ./central/func/${in_root}_smooth.nii.gz -d ./regressors/fsl_des/${bids_root}_regressors_${clean}_des -o ./central/func_${clean}/${bids_root}_b_est --out_res=./central/func_${clean}/${bids_root}_func.nii.gz

# do cleaning with PRF
clean=wPRF
Text2Vest ./regressors/${bids_root}_regressors_${clean}.txt ./regressors/fsl_des/${bids_root}_regressors_${clean}_des
# generate the cleaned data
fsl_glm -i ./central/func/${in_root}_smooth.nii.gz -d ./regressors/fsl_des/${bids_root}_regressors_${clean}_des -o ./central/func_${clean}/${bids_root}_b_est --out_res=./central/func_${clean}/${bids_root}_func.nii.gz

echo "Finished cleaning the data!"
