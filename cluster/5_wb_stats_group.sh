#!/bin/bash

#SBATCH --time=7:55:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=32G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=MYIND_stats_MYPIPEL
#SBATCH --output=out_MYIND-MYPIPEL_MYCONTRAST_group_MYSES-MYTASK.txt
#SBATCH --mail-type=ALL


module load StdEnv/2023  
module load fsl/6.0.7.7 

# choose session and task
ses_num=MYSES
task=MYTASK

# choose autonomic activity index and pipeline
seed=MYIND
pipel=MYPIPEL
contrast=MYCONTRAST

nperm="5000"

# make appropriate directories
cd /home/miedemam/scratch/MGH_derivatives
mkdir ./wb_stats/group/${seed}/
mkdir ./wb_stats/group/${seed}/${pipel}
mkdir ./wb_stats/group/${seed}/${pipel}/${ses_num}
mkdir ./wb_stats/group/${seed}/${pipel}/${ses_num}/${task}
mkdir ./wb_stats/group/${seed}/${pipel}/${ses_num}/${task}/${contrast}

# remove anything in the folder which could mess things up
rm ./wb_stats/group/${seed}/${pipel}/${ses_num}/${task}/${contrast}/*

# now grab all the COPE images, transform, and concatenate
for i in $(seq -f "%03g" 1 22) 
do
sub_num=$(sed -n "${i}p" run_subs_22.txt)
# transform into MNI space
anat_folder=/home/miedemam/scratch/MGH_derivatives/anat/
transform_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-02/rest.feat/
func_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/
if [ "$contrast" = "pos" ]; then
zstat_in=./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}.feat/stats/zstat1.nii.gz
zstat_out=./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}/zstat1_MNI.nii.gz
fi
if [ "$contrast" = "neg" ]; then
zstat_in=./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}.feat/stats/zstat2.nii.gz
zstat_out=./wb_stats/${seed}/${pipel}/${sub_num}/${ses_num}/${task}/zstat2_MNI.nii.gz
fi
ref_1mm=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/CUDA/gcc9/cuda11.0/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz
if [ "$ses_num" = "02" ]; then
applywarp -i $zstat_in -r $ref_1mm -o $zstat_out -w ${transform_root}reg/highres2standard_warp --premat=${func_root}reg/example_func2highres.mat
fi
if [ "$ses_num" = "03" ]; then
# concatenate and apply the necessary transformations
convert_xfm -omat ${func_root}reg/example_func2highres2.mat -concat ${func_root}reg/example_func2highres.mat ${anat_folder}trans_ses03_2_ses02_sub-${sub_num}.mat
applywarp -i $zstat_in -r $ref_1mm -o $zstat_out -w ${transform_root}reg/highres2standard_warp --premat=${func_root}reg/example_func2highres2.mat
fi
# prepare to merge everything
fslmaths $zstat_out -subsamp2 $zstat_out
cp $zstat_out ./wb_stats/group/${seed}/${pipel}/${ses_num}/${task}/${contrast}/sub-${sub_num}_z.nii.gz
done

cd ./wb_stats/group/${seed}/${pipel}/${ses_num}/${task}/${contrast}
fslmerge -t merged_zstats sub*
rm sub*

# now downsample (and smooth)
# fslmaths merged_zstats -subsamp2 merged_zstats_2mm

# now run randomise on merged files
randomise -i merged_zstats -o randomise -1 -x -n $nperm -T
