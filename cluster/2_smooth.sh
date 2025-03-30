#!/bin/bash

#SBATCH --array=1-41
#SBATCH --time=00:28:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=8G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=smooth_over
#SBATCH --mail-type=ALL
#SBATCH --output=out_smooth_%j.txt

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
func_data=${func_root}filtered_func_data.nii.gz
smooth_root=/home/miedemam/scratch/MGH_derivatives/hd_transfer/func/
smooth_data=${smooth_root}sub-${sub_num}_ses-${ses_num}_${task}_smooth
# check if data exists, otherwise, move on:
if [ ! -e "$func_data" ]; then
echo "Functional data not found - moving on."
continue
fi
# apply a 5 mm FWHM smoothing filter with edge preservation
susan $func_data -1 2.123 3 1 0 $smooth_data
# also copy over motion parameters from FEAT
mot_root=/home/miedemam/scratch/MGH_derivatives/hd_transfer/motion/
mot_data=${mot_root}sub-${sub_num}_ses-${ses_num}_${task}_mot.par
cp ${func_root}mc/prefiltered_func_data_mcf.par $mot_data
done
done

echo "Finished preparing data for PRF estimation and transfer."
echo "This took:"
echo $SECONDS
