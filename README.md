# autonomic_preproc
Scripts for processing data for work presented in Miedema, M., Dagenais, R., Torabi, M., Askari, E., Long, S., and Mitsis, G.D. Modelling the effect of cardiac and respiratory fluctuations on the central autonomic network in a novel test-retest dataset. [submitted for publication]

Please note that physiological response function (PRF) estimation and RETROICOR scripts have been adapted from https://github.com/mkassinopoulos/PRF_estimation; for further use please cite Kassinopoulos, M., Mitsis, G.D., 2019. Identification of physiological response functions to correct for fluctuations in resting-state fMRI related to heart rate and respiration. Neuroimage 202, 116150. https://doi.org/10.1016/j.neuroimage.2019.116150

## Guidelines for use related to PRF estimation:
1. Run physio_preproc to prepare PPG and ventilation data for use in PRF estimation, calculating the heart rate (HR) and respiratory flow (RF).
2. Run extract GS_and_CAN to prepare global signal (GS) for use in estimation.
3. Run estimate_group_PRF to estimate PRF parameters for all subjects for a single scan.
4. Run estimate_scan_PRF_constrained to estimate PRF parameters for an individual subject's single scan, starting from the group estimation parameters.
5. If desired, run stability estimation (script 3_PRF_stability in /cluster)

## Guidelines for use related to fMRI data analysis:
1. Perform initial fMRI pre-processing (scripts 1 in /cluster)
2. Run generate_regressors to prepare physiological and motion nuisance regressors for use in denoising
3. Perform subsequent fMRI denoising (script 3_clean_func in /cluster)
4. Run first level statistical analysis (script 5_wb_stats_individ in /cluster)
5. Run group level statistical analysis (script 5_wb_stats_group in /cluster)
