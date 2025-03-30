#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=00:08:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=2G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=mask_over
#SBATCH --mail-type=ALL
#SBATCH --output=out_mask_%j.txt

module load StdEnv/2023  
module load fsl/6.0.7.7 

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs.txt) 

declare -a ses=("02" "03") 
declare -a tasks=("rest" "breathing" "coldpressor")

for ses_num in ${ses[@]}
do
for task in ${tasks[@]}
# no cold pressor task
do
if [ "$ses_num" = "03" ]; then
if [ "$task" = "coldpressor" ]; then
continue
fi
fi
func_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/
# check if data exists, otherwise, move on:
if [ ! -d "$func_root" ]; then
echo "FEAT folder not found - moving on."
continue
fi
# copy over brain mask from FEAT
mask_root=/home/miedemam/scratch/MGH_derivatives/seg/
mask_data=${mask_root}sub-${sub_num}_ses-${ses_num}_${task}_brainmask.nii.gz
cp ${func_root}mask.nii.gz $mask_data
done
done

echo "Finished moving over the brain mask for PRF estimation and transfer."
echo "This took:"
echo $SECONDS
