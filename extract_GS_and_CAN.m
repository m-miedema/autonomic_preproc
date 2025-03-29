%% Extract the GS from fMRI data
% This script extracts 4 variants of GS. Specifically, it extracts the 
% mean timeseries from the whole brain (WB), gray matter (GM),
% white matter (WM) and the cerebrospinal fluid (CSF).
% In addition, this script removes the motion realignment parameters
% as well as the cardiac and respiratory  RETROICOR regressors from 
% each of the 4 variants of GS. In total, 8 variants of GS are extracted
% which are stored in the MAT-File 'GS.mat'.

%% 1: Define directories for necessary data

clear, clc, close all


%  Set the following parameters !!
% 
% sub_num = '0059';
% ses_num = "03";
% task = "breathing"; 

subjects = {'0059','0704','0134','1092','1687','1693','1750','3046','4901','5262','5571','5931',...
    '6293','7002','7712','7994','8116','8223','8420','8803','8907','9922'};

for sub_i = 1:length(subjects)

    sub_num=subjects{sub_i};

    % for 7994, there's an inconsistent number of triggers

    for i = 1:3
        
        % convenient designation of sessions and scans
        if i == 1
            ses_num="02";
            task="rest";
        elseif i == 2
            ses_num="02";
            task="breathing";
        elseif i == 3
            ses_num="02";
            task="coldpressor";
        elseif i == 4
            ses_num="03";
            task="rest";
        elseif i == 5
            ses_num="03";
            task="breathing";
        elseif i == 6
            ses_num="03";
            task="visual";
        end
        
        fprintf('    Now processing sub-%s, ses-%s, %s. ',sub_num,ses_num,task)
        
        % -----------------------------------------
        
        % define path to root directory
        root_dir = "D:\estimation";
        scan_root = strcat("sub-",sub_num,"_ses-",ses_num,"_",task);
        
        path_func = strcat(root_dir,'\func\',scan_root,"_smooth.nii.gz");
        path_mask = strcat(root_dir,'\seg\',scan_root,'_brainmask.nii.gz');
        path_CSF = strcat(root_dir,'\seg\',scan_root,'_CSF.nii.gz');
        path_GM = strcat(root_dir,'\seg\',scan_root,'_GM.nii.gz');
        path_WM = strcat(root_dir,'\seg\',scan_root,'_WM.nii.gz');
        path_CAN = strcat(root_dir,'\seg\',scan_root,'_CAN.nii.gz');
        path_movement = strcat(root_dir,'\motion\',scan_root,'_mot.par');
        
        bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/sub-",sub_num,"_ses-",ses_num,"_task-",task);
        physio_out_dir = "D:\estimation\physio_preproc\";
        path_physio = strcat(physio_out_dir,bids_root,'_physio_and_triggers.mat');
        
        out_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
        path_output = strcat(root_dir,'\extracted\',out_root,'_GS.mat');
        path_CAN_output = strcat(root_dir,'\extracted\',out_root,'_CAN.mat');
        
        %%  2: Load niftii files
        
        % nii_func = load_nii(path_func); img = nii_func.img;
        % nii_mask = load_nii(path_mask);   img_mask = nii_mask.img;
        % nii_CSF = load_nii(path_CSF);   img_CSF = nii_CSF.img;
        % nii_GM = load_nii(path_GM);   img_GM = nii_GM.img;
        % nii_WM = load_nii(path_WM);  img_WM = nii_WM.img;
        % nii_CAN = load_nii(path_CAN);  img_CAN = nii_CAN.img;
        img = double(niftiread(path_func));
        img_mask = double(niftiread(path_mask));
        img_CSF = double(niftiread(path_CSF));
        img_GM = double(niftiread(path_GM));
        img_WM = double(niftiread(path_WM));
        img_CAN = double(niftiread(path_CAN));
        [NX, NY, NZ, NV] = size(img);
        
        %% 3: Copy all voxel timeseries from the whole-brain in a 2-d matrix
        
        NVX = 0;
        CordAll=zeros(NX*NY*NZ,3);
        for x=1:NX
            for y=1:NY
                for z=1:NZ
                    if img_mask(x,y,z)==1
                        voxel = squeeze(img(x,y,z,:));
                        voxel_mean = mean(voxel); 
                        if   voxel_mean > 100
                            NVX=NVX+1;
                            CordAll(NVX,1)=x;
                            CordAll(NVX,2)=y;
                            CordAll(NVX,3)=z;
                            img(x,y,z,:) = detrend(voxel,'linear') + voxel_mean; 
                        end                
                    end
                end
            end
        end
        CordAll=CordAll(1:NVX,:);
        
        img_col=zeros(NV,NVX);
        for vox=1:NVX
            x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3); 
            img_col(:,vox) = squeeze(img(x,y,z,:));
        end
        
        %% 4: Extract the mean timeseries from WB, GM, WM and CSF
        
        idx_GM=zeros(NVX,1); N_GM=0;
        idx_WM=zeros(NVX,1); N_WM=0;
        idx_CSF=zeros(NVX,1); N_CSF=0;
        for vox=1:NVX
            x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
            if img_GM(x,y,z)>0.5
                N_GM=N_GM+1;
                idx_GM(N_GM)=vox;
            elseif img_CSF(x,y,z)>0.8
                N_CSF=N_CSF+1;
                idx_CSF(N_CSF)=vox;
            elseif img_WM(x,y,z)>0.9
                N_WM=N_WM+1;
                idx_WM(N_WM)=vox;
            end
        end
        
        idx_GM=idx_GM(1:N_GM,1);
        idx_WM=idx_WM(1:N_WM,1);
        idx_CSF=idx_CSF(1:N_CSF,1);
        fprintf('   GM:  %3.1f%%  \n',100*N_GM/NVX);
        fprintf('   WM:  %3.1f%%  \n',100*N_WM/NVX);
        fprintf('   CSF: %3.1f%%  \n',100*N_CSF/NVX);
        
        CSF = mean(img_col(:,idx_CSF)')';
        GM = mean(img_col(:,idx_GM)')';
        WM = mean(img_col(:,idx_WM)')';
        WB = mean(img_col')';
        
        fprintf(' ----- Correlation ----- \n')
        fprintf('CSF-GM: %3.1f%%  \n',100*corr(CSF,GM))
        fprintf('WM-GM: %3.1f%%  \n',100*corr(WM,GM))
        fprintf('CSF-WM: %3.1f%%  \n',100*corr(CSF,WM))
        fprintf('WB-GM: %3.1f%%  \n',100*corr(WB,GM))
        
        TR = 1.03;   timeMR = 0 :TR: (NV-1)*TR;
        % 
        % plot(timeMR,[(GM),(WM),(CSF)])
        % xlabel('Time (s)')
        % legend('Gray matter','White matter','Cerebrospinal fluid')
        
        %% 5:  Create regressor matrix with motion parameters and RETROICOR regressors
        
        movRegr = zscore(load(path_movement));
        mov_d_Regr = diff(movRegr);
        mov_d_Regr = [mov_d_Regr(1,:);mov_d_Regr];
        load(path_physio);
        
        ind_BOLD = find(trig_ind==1);
        addpath('RETROICOR')
        RETR_RespRegr = func_RETR_Resp_regressors(resp,3,Fs,rsp_phase_interp); NR_RETRresp=size(RETR_RespRegr,2); RETR_RespRegr = RETR_RespRegr(ind_BOLD,:);
        RETR_CardRegr = func_RETR_Card_regressors(time,PPGlocs,3);   NR_card=size(RETR_CardRegr,2);  RETR_CardRegr = RETR_CardRegr(ind_BOLD,:);
        
        % drop first volumes
        dropN = 3;
        regr = [RETR_RespRegr(dropN+1:end,:), RETR_CardRegr(dropN+1:end,:)];    
        regr =  detrend(regr,'linear');
        rLinear = [1:NV]';
        regr = [regr, movRegr, mov_d_Regr, rLinear, ones(NV,1)];
        
        %% 6: Remove confounds from the 4 GS variants extracted in Section 4 and save variables
        
        B = regr\GM;  GM_clean = GM - regr*B;
        B = regr\WM;  WM_clean = WM - regr*B;
        B = regr\CSF;  CSF_clean = CSF - regr*B;
        B = regr\WB;  WB_clean = WB - regr*B;
        
        %% 7: Extract mean signal from ROIs associated with the CAN
        CAN_mean = zeros(NV,19);
        CAN_clean_mean = zeros(NV,19);
        for k=1:19
            idx_CAN=zeros(NVX,1); N_CAN=0;
            for vox=1:NVX
                x=CordAll(vox,1);    y=CordAll(vox,2);    z=CordAll(vox,3);
                if img_CAN(x,y,z,k)>0.5
                    N_CAN=N_CAN+1;
                    idx_CAN(N_CAN)=vox;
                end
            end
            
            idx_CAN=idx_CAN(1:N_CAN,1);    
            CAN_mean(:,k) = mean(img_col(:,idx_CAN)')';
            CAN_clean_mean(:,k) = CAN_mean(:,k) - regr*B;
        
        end
        
        %% 8: Export results
        save(path_output,'WB','GM','WM','CSF','WB_clean','GM_clean','WM_clean','CSF_clean')
        save(path_CAN_output,'CAN_mean','CAN_clean_mean')

        clearvars -except subjects sub_i sub_num i
    end
end






























