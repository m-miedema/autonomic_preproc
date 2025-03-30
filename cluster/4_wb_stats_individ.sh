#!/bin/bash

#SBATCH --array=1-22
#SBATCH --time=0:38:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=16G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=MYIND_stats_MYPIPEL
#SBATCH --mail-type=ALL
#SBATCH --output=out_MYIND_individ_MYSES-MYTASK_%j.txt

module load StdEnv/2023  
module load fsl/6.0.7.7 

# choose session and task
ses_num=MYSES
task=MYTASK
nVol=MYVOL

# choose autonomic activity index and pipeline
seed=MYIND
pipel=MYPIPEL

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_22.txt) 

echo "Processing subject:"
bids_root=sub-${sub_num}_ses-${ses_num}_task-${task}
echo $bids_root

echo "for pipeline:"
echo $pipel

# truncate scan by 20 volumes
# fslroi ./central/func_${pipel}/${bids_root}_func ./central/func_${pipel}/${bids_root}_func_dropped 20 $nVol

# make appropriate directories

mkdir ./wb_stats/${seed}/
mkdir ./wb_stats/${seed}/${pipel}
mkdir ./wb_stats/${seed}/${pipel}/${sub_num}
mkdir ./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}
mkdir ./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}

cp ./wb_1lev_mod_no_reg.fsf ./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}/stats1.fsf

cd ./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}/

sed -i "s|subXXX|${sub_num}|g" stats1.fsf
sed -i "s|sesXXX|${ses_num}|g" stats1.fsf
sed -i "s|taskXXX|${task}|g" stats1.fsf
sed -i "s|pipeXXX|${pipel}|g" stats1.fsf
sed -i "s|nuXXX|${seed}|g" stats1.fsf
sed -i "s|volXXX|${nVol}|g" stats1.fsf

feat stats1.fsf
