%% Generate connectivity edges within the CAN

clear, clc, close all


subjects = {'0059','0704','0134','1092','1687','1693','1750','3046','4901','5262','5571','5931',...
    '6293','7002','7712','7994','8116','8223','8420','8803','8907','9922'};

tasks = {'rest','breathing','cold-pressor','rest','breathing'};
sessions = {'2','2','2','3','3'};
c_methods = {'minimal' 'wPRF'};

% define paths
root_dir = "D:\estimation";

N_s = length(subjects);

for c_type = 1:2
    cleaning = c_methods{c_type};
    all_CAN_conn = {};
    mean_CAN_conn = {};

    for i = 1:5
    
        CAN_conn_i = zeros(19,19,N_s);
        
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
        end
    
        for sub_i = 1:N_s
    
            sub_num=subjects{sub_i};
            
            fprintf('    Now processing sub-%s, ses-%s, %s. \n',sub_num,ses_num,task)
            
            % ----------------------------------------- LOAD DATA
    
            scan_root = strcat("sub-",sub_num,"_ses-",ses_num,"_",task);
            bids_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
    
            path_CAN = strcat(root_dir,'\seg\',scan_root,'_CAN.nii.gz');
            path_func = strcat(root_dir,'\func_',cleaning,'\',bids_root,'_func.nii.gz');
            path_conn_output = strcat(root_dir,'\CAN_conn\all\',bids_root,'_',cleaning,'_CAN_conn_z.mat');
    
            img = double(niftiread(path_func));
            img_CAN = double(niftiread(path_CAN));
            [NX, NY, NZ, NV] = size(img);
    
            % extract mean signal in CAN ROIs
            CAN_mean = zeros(NV,size(img_CAN,4));
            for j = 1:size(img_CAN,4)
                ROI_j = repmat(img_CAN(:,:,:,j)>0,[1 1 1 NV]);
                mean_j = reshape(mean(ROI_j.*img,[1 2 3]),1,[]);
                CAN_mean(:,j) = mean_j;
            end
    
            % calculate connectivity and z-score
            sCAN_conn = corrcoef(CAN_mean);
            CAN_conn_i(:,:,sub_i) = atanh(sCAN_conn);
    
            if sum(find(isnan(CAN_conn_i))) > 0
                disp("Caution, CAN signals missing: check mask definition.")
            end
    
            save(path_conn_output,"CAN_mean","sCAN_conn")
        end
    
        all_CAN_conn{end+1} = CAN_conn_i;
        mean_CAN_conn{end+1} = mean(abs(CAN_conn_i),3);
    end
    
    path_all_conn_output = strcat(root_dir,'\CAN_conn\',cleaning,'_CAN_conn_z.mat');
    save(path_all_conn_output,"all_CAN_conn")
end

% create some figures
cleaning='wPRF';
path_all_conn_output = strcat(root_dir,'\CAN_conn\',cleaning,'_CAN_conn_z.mat');
load(path_all_conn_output,"all_CAN_conn")

for i = 1:5

    CAN_conn_i = mean(abs(all_CAN_conn{i}),3);
    
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
    end

    out_root = strcat(cleaning,"_ses-",ses_num,"_task-",task);
    conn_fig_path = strcat(root_dir,'\CAN_conn\',out_root,'_mean_CAN_conn_z.png');


    % prep for correlation plot
    [colormap2]=cbrewer2('div','RdBu',100);
    %close all
    % now plot!
    fig=figure(2);
    imshow(CAN_conn_i,[0 1],'InitialMagnification','fit')%[-0.2 0.2],
    hold on
    axis on
    set(gca,'TickLength',[0 0])
    xlabel('ROIs in the CAN')
    ax = gca;
    ax.FontSize = 7; 
    colormap(flipud(colormap2))
    cb = colorbar;
    cb.Label.String = 'Mean absolute correlation';
    cb.Label.FontAngle = 'italic';
    title(sprintf('Mean CAN connectivity for %s (session-%s) ',tasks{i},sessions{i}),'Interpreter','none')
    saveas(fig,conn_fig_path)
end

































