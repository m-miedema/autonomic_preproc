#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=0:48:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=16G  # memory per CPU core  
#SBATCH --job-name=transforms
#SBATCH --mail-type=ALL
#SBATCH --output=out_trans_%j.txt

module load StdEnv/2023  
module load fsl/6.0.7.7 

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_22.txt) 

echo "Processing subject:"
echo $sub_num

# first, calculate transforms between T1s images from separate sessions
anat_folder=/home/miedemam/scratch/MGH_derivatives/anat/
anat_file2=${anat_folder}sub-${sub_num}_ses-02_acq-mprage_T1w_brain
anat_file3=${anat_folder}sub-${sub_num}_ses-03_acq-mprage_T1w_brain
flirt -in $anat_file2 -ref $anat_file3 -omat ${anat_folder}trans_ses02_2_ses03_sub-${sub_num}.mat
flirt -in $anat_file3 -ref $anat_file2 -omat ${anat_folder}trans_ses03_2_ses02_sub-${sub_num}.mat

# next, transform anatomical segmentations for each session and scan

declare -a ses=("02" "03") 
declare -a tasks=("rest" "breathing" "coldpressor")

for ses_num in ${ses[@]}
do
for task in ${tasks[@]}
do
echo "Processing session:"
echo $ses_num
echo "Processing task:"
echo $task
# segmentation was performed on session 2 data only
seg_root_in=${anat_folder}sub-${sub_num}_ses-02_acq-mprage_T1w_brain_pve
seg_root_out=/home/miedemam/scratch/MGH_derivatives/seg/sub-${sub_num}_ses-${ses_num}_${task}
func_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/
func_data=${func_root}filtered_func_data.nii.gz
# check if data exists, otherwise, move on:
if [ ! -e "$func_data" ]; then
echo "Functional data not found - moving on."
continue
fi
if [ "$ses_num" = "02" ]; then
flirt -in ${seg_root_in}_0.nii.gz -ref $func_data -applyxfm -init ${func_root}reg/example_func2highres_inv.mat -out ${seg_root_out}_CSF.nii.gz
flirt -in ${seg_root_in}_1.nii.gz -ref $func_data -applyxfm -init ${func_root}reg/example_func2highres_inv.mat -out ${seg_root_out}_GM.nii.gz
flirt -in ${seg_root_in}_2.nii.gz -ref $func_data -applyxfm -init ${func_root}reg/example_func2highres_inv.mat -out ${seg_root_out}_WM.nii.gz
fi
if [ "$ses_num" = "03" ]; then
# concatenate the necessary transformations
convert_xfm -omat ${func_root}reg/example_func2highres_inv2.mat -concat ${func_root}reg/example_func2highres_inv.mat ${anat_folder}trans_ses02_2_ses03_sub-${sub_num}.mat
flirt -in ${seg_root_in}_0.nii.gz -ref $func_data -applyxfm -init ${func_root}reg/example_func2highres_inv2.mat -out ${seg_root_out}_CSF.nii.gz
flirt -in ${seg_root_in}_1.nii.gz -ref $func_data -applyxfm -init ${func_root}reg/example_func2highres_inv2.mat -out ${seg_root_out}_GM.nii.gz
flirt -in ${seg_root_in}_2.nii.gz -ref $func_data -applyxfm -init ${func_root}reg/example_func2highres_inv2.mat -out ${seg_root_out}_WM.nii.gz
fi
echo "Finished transforming segmentations!"
done
done

# note that these masks are NOT binarized but remain probabalistic
echo "Finished all transformations to functional space!"
echo "This took:"
echo $SECONDS
