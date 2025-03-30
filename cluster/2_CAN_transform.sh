#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=0:58:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=16G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=CAN_mask_transform
#SBATCH --mail-type=ALL
#SBATCH --output=out_CAN_mask_%j.txt

module load StdEnv/2023  
module load fsl/6.0.7.7 

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs.txt) 

echo "Processing subject:"
echo $sub_num

# original CAN mask in MNI152 space
CAN_MNI1mm=/home/miedemam/scratch/MGH_derivatives/CAN_mask/CAN_masks.nii.gz

# only session 2, rest scan has nonlinear transformation to MNI152 space
transform_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-02/rest.feat/

# calculate an inverted transformation from MNI152 to session 2 structural
#invwarp -w ${transform_root}reg/highres2standard_warp -o ${transform_root}reg/standard2highres_warp -r ${transform_root}reg/highres


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
seg_out=/home/miedemam/scratch/MGH_derivatives/seg/sub-${sub_num}_ses-${ses_num}_${task}_CAN
func_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/
# check if data exists, otherwise, move on:
if [ ! -d "$func_root" ]; then
echo "FEAT directory not found for this scan - moving on."
continue
fi
if [ "$ses_num" = "02" ]; then
applywarp -i $CAN_MNI1mm -r ${func_root}filtered_func_data.nii.gz  -o $seg_out -w ${transform_root}reg/standard2highres_warp --postmat=${func_root}reg/example_func2highres_inv.mat
fi
if [ "$ses_num" = "03" ]; then
# concatenate the necessary transformations
applywarp -i $CAN_MNI1mm -r ${func_root}filtered_func_data.nii.gz  -o $seg_out -w ${transform_root}reg/standard2highres_warp --postmat=${func_root}reg/example_func2highres_inv2.mat
fi
echo "Finished transforming CAN mask!"
# now binarize the mask
fslmaths $seg_out -thr 0.05 -bin $seg_out
done
done

echo "Finished creating all masks in functional space!"
echo "This took:"
echo $SECONDS
