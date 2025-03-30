#!/bin/bash

#SBATCH --array=1-5
#SBATCH --time=02:58:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=16G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=FEAT_cold_new
#SBATCH --mail-type=ALL

module load StdEnv/2023  
module load fsl/6.0.7.7 

# choose the task to process
task=coldpressor
ntp=780

# choose the session
session=02

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_Feb.txt) 

sub_folder=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/
ses_folder=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${session}/

template_file=/home/miedemam/scratch/MGH_derivatives/feat/preproc_template.fsf


# note that if already run, no need to repeat nonlinear registration
# ideally, adjust job time accordingly before submission
full_reg=0

if [ ! -d "$ses_folder" ]; then
  echo "${ses_folder} does not exist, therefore FEAT has not been run for this subject/session. The folder will be created and nonlinear registration will be performed."
  mkdir $sub_folder
  mkdir $ses_folder
  full_reg=1
  
  # also need to copy over the non-BET structural image to keep FEAT happy in this case
  anat_folder=/home/miedemam/scratch/bids_MGH_FINAL/sub-${sub_num}/ses-${session}/anat/
  bet_folder=/home/miedemam/scratch/MGH_derivatives/anat/
  cp ${anat_folder}sub-${sub_num}_ses-${session}_acq-mprage_T1w.nii.gz ${bet_folder}sub-${sub_num}_ses-${session}_acq-mprage_T1w.nii.gz 
fi

echo "Processing subject:"
echo $sub_num
echo "Processing subject:"
echo $session
echo "Processing task:"
echo $task

# modify the template file
cp $template_file ${ses_folder}/preproc_${task}.fsf
cd $ses_folder
sed -i "s|subXXX|${sub_num}|g" preproc_${task}.fsf
sed -i "s|sesXXX|${session}|g" preproc_${task}.fsf
sed -i "s|taskXXX|${task}|g" preproc_${task}.fsf
sed -i "s|tpXXX|${ntp}|g" preproc_${task}.fsf
sed -i "s|regmniXXX|${full_reg}|g" preproc_${task}.fsf

# run FEAT
feat preproc_${task}.fsf

echo "Finished basic preprocessing with FEAT!"
echo "This took:"
echo $SECONDS
