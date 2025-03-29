% compares outcomes for PRF estimation across numerous iterations

clear, close all, clc
%% set up loop across subjects
subjects = {'0059','0704','0134','1092','1687','1693','1750','3046','4901','5262','5571','5931',...
    '6293','7002','7712','7994','8116','8223','8420','8803','8907','9922'};
% input directory
root_dir = "D:\estimation\PRF_stability\";

Ts_10=0.1;
t_IR = 0:Ts_10:(600*Ts_10);
N_iter = 25;

load('condition_names.mat')

CRF_stds = zeros(length(subjects),5,length(t_IR));
RRF_stds = zeros(length(subjects),5,length(t_IR));

% assess
yPred_sim_r = zeros(length(subjects),10);
yPred_sim_c = zeros(length(subjects),10);

for sub_i = 1:length(subjects)

    sub_num=subjects{sub_i};
    
    CRF_all = {};%zeros(N_iter,601,5);
    RRF_all = {};%zeros(N_iter,601,5);
    yPred_all = {};%zeros(N_iter,NV,5);
    yPred_all_r = {};
    yPred_all_c = {};

    HR_all = {};
    RF_all = {};
    trig_ind_all = {};

    for i = 1:5
        
        % convenient designation of sessions and scans
        if i == 1
            ses_num="02";
            task="rest";
            NV = 297;
        elseif i == 2
            ses_num="02";
            task="breathing";
            NV = 397;
        elseif i == 3
            ses_num="02";
            task="coldpressor";
            NV = 777;
        elseif i == 4
            ses_num="03";
            task="rest";
            NV = 597;
        elseif i == 5
            ses_num="03";
            task="breathing";
            NV = 397;
        end

        yPred = zeros(N_iter,NV);
        yPred_r = zeros(N_iter,NV);
        yPred_c = zeros(N_iter,NV);
        CRF_all_i = zeros(N_iter,601);
        RRF_all_i = zeros(N_iter,601);
        
        fprintf('    Now processing sub-%s, ses-%s, %s. \n',sub_num,ses_num,task)

        % load in HR, RF                
        ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
        bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/",ext_root);
        physio_out_dir = "D:\estimation\physio_preproc\";
        path_physio = strcat(physio_out_dir,bids_root,'_physio_and_triggers.mat');
            
        load(path_physio,'HR_10','RF','trig_ind_10')
        HR_all{end+1} = HR_10;
        RF_all{end+1} = RF;
        trig_ind_all{end+1} = trig_ind_10;

        for iter = 1:N_iter
            % load in the data for this iteration
            ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
            path_out = strcat(root_dir,'/iter_',num2str(iter,'%03.f'),'/',ext_root,'_PRF.mat');
            load(path_out,'CRF_sc','RRF_sc','r_PRF_sc','yPred','yPred_card','yPred_resp','CRF_par','RRF_par','scale_par','xPhys')

            % store in new arrays
            CRF_all_i(iter,:) = CRF_sc;
            RRF_all_i(iter,:) = RRF_sc;

            % might need to check/flip the sign
            B=xPhys\yPred;
            if B(1) < 0
                CRF_all_i(iter,:) = -CRF_sc;
            end
            if CRF_all_i(iter,20) < 0
                CRF_all_i(iter,:) = -CRF_all_i(iter,:);
            end
            if B(2) < 0
                RRF_all_i(iter,:) = -RRF_sc;
            end
        
        end

        CRF_all{end+1} = CRF_all_i;
        RRF_all{end+1} = RRF_all_i;
        yPred_all{end+1} = yPred;
        yPred_all_r{end+1} = yPred_r;
        yPred_all_c{end+1} = yPred_c;
    end

    % now assess inter-scan PRF agreement
    
    for j = 1:5
        if j == 1
            ind1 = 1;
            ind2 = 4;
            j1 = 1;
            j2 = 2;
        elseif j == 2
            ind1 = 2;
            ind2 = 5;
            j1 = 3;
            j2 = 4;
        elseif j == 3
            ind1 = 1;
            ind2 = 3;
            j1 = 5;
            j2 = 6;
        elseif j == 4
            ind1 = 1;
            ind2 = 2;
            j1 = 7;
            j2 = 8;
        elseif j == 5
            ind1 = 4;
            ind2 = 5;
            j1 = 9;
            j2 = 10;
        end
        % yPred1 = yPred_all{ind1};
        % yPred2 = yPred_all{ind2};
        % yPred1_r = yPred_all_r{ind1};
        % yPred2_r = yPred_all_r{ind2};
        % yPred1_c = yPred_all_c{ind1};
        % yPred2_c = yPred_all_r{ind2};

        HR1 = HR_all{ind1};
        HR2 = HR_all{ind2};
        RF1 = RF_all{ind1};
        RF2 = RF_all{ind2};
        trig_ind1 = find(trig_ind_all{ind1});
        trig_ind2 = find(trig_ind_all{ind2});

        CRF1 = CRF_all{ind1};
        CRF2 = CRF_all{ind2};
        RRF1 = RRF_all{ind1};
        RRF2 = RRF_all{ind2};

        corr12_c = zeros(N_iter,1);
        corr12_r = zeros(N_iter,1);
        corr21_c = zeros(N_iter,1);
        corr21_r = zeros(N_iter,1);

        % estimate sessions with opposing and original PRFs
        for iter = 1:N_iter
            CRF1_iter = CRF1(iter,:);
            RRF1_iter = RRF1(iter,:);
            CRF2_iter = CRF2(iter,:);
            RRF2_iter = RRF2(iter,:);

            HR_conv12 = conv(HR1,CRF2_iter); HR_conv12_MR = HR_conv12(trig_ind1);
            RF_conv12 = conv(RF1,RRF2_iter); RF_conv12_MR = RF_conv12(trig_ind1);

            HR_conv21 = conv(HR2,CRF1_iter); HR_conv21_MR = HR_conv21(trig_ind2);
            RF_conv21 = conv(RF2,RRF1_iter); RF_conv21_MR = RF_conv21(trig_ind2);

            HR_conv11 = conv(HR1,CRF1_iter); HR_conv11_MR = HR_conv11(trig_ind1);
            RF_conv11 = conv(RF1,RRF1_iter); RF_conv11_MR = RF_conv11(trig_ind1);

            HR_conv22 = conv(HR2,CRF2_iter); HR_conv22_MR = HR_conv22(trig_ind2);
            RF_conv22 = conv(RF2,RRF2_iter); RF_conv22_MR = RF_conv22(trig_ind2);

            % calculate correlations
            corr12_c(iter) = corr(HR_conv12_MR.',HR_conv11_MR.','rows','complete');
            corr12_r(iter) = corr(RF_conv12_MR,RF_conv11_MR,'rows','complete');
            corr21_c(iter) = corr(HR_conv21_MR.',HR_conv22_MR.','rows','complete');
            corr21_r(iter) = corr(RF_conv21_MR,RF_conv22_MR,'rows','complete');
        end

        % save mean correlations
        yPred_sim_c(sub_i,j1) = mean(abs(corr12_c),'omitnan');
        yPred_sim_c(sub_i,j2) = mean(abs(corr21_c),'omitnan');
        yPred_sim_r(sub_i,j1) = mean(abs(corr12_r),'omitnan');
        yPred_sim_r(sub_i,j2) = mean(abs(corr21_r),'omitnan');
    end


end

sim_c = reshape(yPred_sim_c,[2*length(subjects) 5]);
sim_r = reshape(yPred_sim_r,[2*length(subjects) 5]);
save("intrascan_PRF_prediction_similarity.mat","sim_r","sim_c")

summaryfig = figure;
t1 = tiledlayout(1,1,'TileSpacing','compact');
group_names = {"Rest 1 x 2","Breathing 1 x 2", "Rest 1 x Cold 1","Rest 1 x Breathing 1","Rest 2 x Breathing 2"};
comp_names = {"Cardiac" "Respiration"};
nexttile
title('Cardiac')
h = daboxplot({sim_c,sim_r},'scatter',2,'whiskers',0,'boxalpha',0.7,...
            'xtlabels', group_names,'legend',comp_names); 
h.lg.Layout.Tile = 'south';
ylabel('Correlation')
title('PRF generalizability across scans')

% 
% summaryfig = figure;
% t1 = tiledlayout(1,1,'TileSpacing','compact');
% group_names = {"Rest 1 x 2","Breathing 1 x 2", "Rest 1 x Cold 1","Rest 1 x Breathing 1","Rest 2 x Breathing 2"};
% comp_names = {"Cardiac" "Respiration"};
% nexttile
% title('Cardiac')
% h = daboxplot({yPred_sim_c(:,1:2:end),yPred_sim_c(:,2:2:end),yPred_sim_r(:,1:2:end),yPred_sim_r(:,2:2:end)},'scatter',2,'whiskers',0,'boxalpha',0.7,...
%             'xtlabels', group_names,'legend',comp_names); 
% h.lg.Layout.Tile = 'south';
% ylabel('Correlation')
% title('PRF generalizability across scans')