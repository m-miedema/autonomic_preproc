#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=0:39:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=4G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=stability-MYITER
#SBATCH --mail-type=ALL
#SBATCH --output=out_stability-MYITER_%j.txt

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_22.txt) 

# set the iteration
iteration=MYITER

echo "Processing subject:"
echo $sub_num

module load matlab/2023b.2

echo "Matlab loaded! Running estimation script:"

matlab -nodisplay -r "estimate_scan_PRF_stability(${sub_num},${iteration})" -sd /home/miedemam/scratch/MGH_derivatives/matlab_functions

echo "Finished estimating PRFs for all scanning segments!"
echo "This took:"
echo $SECONDS
