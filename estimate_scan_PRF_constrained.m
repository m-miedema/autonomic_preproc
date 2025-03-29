
clear, close all, clc
%% set up loop across subjects
subjects = {'0059','0704','0134','1092','1687','1693','1750','3046','4901','5262','5571','5931',...
    '6293','7002','7712','7994','8116','8223','8420','8803','8907','9922'};

for sub_i = 1:length(subjects)

    sub_num=subjects{sub_i};

    % for 7994, there's an inconsistent number of triggers - robust to
    % selecting final scans

    % parameters I could potentially set
    % InitialPopulationRange

    for i = 1:5
        
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
        
        fprintf('    Now processing sub-%s, ses-%s, %s. \n',sub_num,ses_num,task)
        


        %%  1:  Load physiological variables (heart rate and respiration) and global signal (GS) from MAT-File
                
        % define path to root directory
        root_dir = "D:\estimation";
        ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
        bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/",ext_root);
        
        physio_out_dir = "D:\estimation\physio_preproc\";
        path_physio = strcat(physio_out_dir,bids_root,'_physio_and_triggers.mat');
        path_GS = strcat(root_dir,'\extracted\',ext_root,'_GS.mat');
        
        path_out = strcat(root_dir,'\PRF_constrained_nofilt\',ext_root,'_PRF.mat');
        
        % -----------------------------------------
        
        load(path_physio,"HR_10","resp_10","RF","trig_ind_10","time_10")
        load(path_GS,"WB_clean")
        GS = WB_clean; % test this
        
        Ts_10 = 0.1 ;                                                       % Sampling period in seconds
        %time_10 = 0:Ts_10:(length(HR)-1)*Ts_10;
       

        NV = length(GS);
        % dropping first few volumes
        timeMR = time_10(find(trig_ind_10,NV,"last"));
        trig_ind_10 = find(trig_ind_10,NV,"last");
        % deal with NaNs until first heart beat happens
        firstHR = find(~isnan(HR_10),1,"first");
        lastHR = find(~isnan(HR_10),1,"last");
        HR_10(1:firstHR)=HR_10(firstHR);
        HR_10(lastHR:end)=HR_10(lastHR);
        HR_10 = HR_10.';
        
        % figure()
        % ax1 = subplot(3,1,1);
        % plot(time_10,HR_10)
        % title('Heart rate (HR)')
        % ylabel('HR (bpm)')
        % 
        % ax2 = subplot(3,1,2);
        % plot(time_10,resp_10)
        % title('Respiration')
        % ylabel('Amplitude (a.u.)')
        % 
        % ax3 = subplot(3,1,3);
        % plot(timeMR,GS);
        % title('Global signal (GS)')
        % ylabel('Amplitude (a.u.)')
        % xlabel('Time (s)')
        % linkaxes([ax1,ax2,ax3],'x')
        % xlim([0,max(time_10)])
        
        %% 2: Estimate PRF parameters
        
        %resp_s = smooth(resp,10*1.5) ;
        %RF=diff(resp_s); RF=[0;RF(:)]; RF = RF.^2;

        filt_all = false;
        if filt_all
        
            fs_MR = 1/1.03;
            fs_phys = 10;
            lp_filt = 0.15;
            [filt_MR1, filt_MR2] = butter(2,lp_filt/(fs_MR/2));
            GS = filtfilt(filt_MR1,filt_MR2,GS);
            [filt_phys1, filt_phys2] = butter(2,lp_filt/(fs_phys/2));
            HR_10 = filtfilt(filt_phys1,filt_phys2,HR_10);
            RF = filtfilt(filt_phys1,filt_phys2,RF);
        end

        
        ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',400,'Display','final','PopulationSize',400);%'UseParallel',1);   % Display: iter
        options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
            'UseParallel',false,'MaxIterations',400,'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-9);%,'PlotFcn','optimplotfval');    % 'PlotFcn','optimplotfval'
        
        %PRF_par = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5 ];
        path_group = strcat(task,'_',ses_num,'_nofilt.mat');
        load(path_group,"PRF_par")
        %PRF_par = PRF_par(1:8);
        PRF_par_group = PRF_par;
        ub = PRF_par+3;
        lb = PRF_par-3; lb(find(lb<0))=0.0001;
        lb(9:10) = -inf; ub(9:10)=inf;
        % lb(5) = 0.8;
        % lb(6) = 0.2;
        % lb(7) = 0.8;
        % lb(8) = 0.2;
        %ub = PRF_par+10;
        %lb = PRF_par-6; lb(find(lb<0))=0;
        %plotobjective(@shufcn,[ub; lb]);

        h_train = @(P) func_M4_PRF_sc(P,Ts_10,HR_10,RF,trig_ind_10,GS,1);
        h_constraint = @bound_constraint;

        % Uncomment the following line if you want to use  Genetic Algorithm
        % (GA). GA may yield better fit with the cost of longer computational time.
        PRF_par = ga(h_train,length(ub),[],[],[],[],lb,ub,h_constraint,[],ga_opts);
        %PRF_par = fmincon(h_train,PRF_par,[],[],[],[],lb,ub,[],options);
        
        h = @(P) func_M4_PRF_sc(P,Ts_10,HR_10,RF,trig_ind_10,GS,0);
        [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred,yPred_card,yPred_resp, HR_conv_MR, RF_conv_MR] = h(PRF_par);
        
        %fprintf('    The sub-%s, ses-%s, %s. \n',sub_num,ses_num,task)

        fprintf(' ----------------------------------------------- \n')
        fprintf('Correlation b/w GS and PRF output \n')
        fprintf('CRF (HR): %3.2f  \n',r_PRF_sc(2))
        fprintf('RRF (RF): %3.2f  \n',r_PRF_sc(3))
        fprintf('CRF & RRF (HR & RF): %3.2f  \n',r_PRF_sc(1))
        
        %%  3: Plot output of PRF model (timeseries and PRF curves)  
        
        %  Set the following parameters !!
        
        smoothPar = 5;
        fontTitle = 20;
        fontLabels = 8;
        fontTxt = 16;
        lineW = 3;
        yl1 = -5.3; yl2 = 5.5;
        
        % -----------------------------------------
        
        t_IR = 0:Ts_10:(length(CRF_sc)-1)*Ts_10;
        
        % screenSize = get(0,'ScreenSize'); xL = screenSize(3); yL = screenSize(4);
        % figure
        % set(gcf, 'Position', [0.2*xL 0.2*yL  0.6*xL 0.6*yL ]);
        % set(gcf, 'Position', [0.1*xL 0.1*yL  0.8*xL 0.8*yL ]);
        % 
        % ax1 = subplot(5,3,1:2);
        % plot(time_10,HR_10)
        % ylabel('HR (bpm)')
        % title(sprintf('Heart rate (HR; %2.0f±%1.0f bpm )',mean(HR_10),std(HR_10)))
        % 
        % ax6 = subplot(5,3,[3,6]);
        % plot(t_IR,CRF_sc,'LineWidth',4), grid on
        % title('Cardiac Response Function (CRF_{sc}) ')
        % xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
        % xlim([0 60])
        % 
        % ax2 = subplot(5,3,4:5);
        % h1=plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
        % h2=plot(time_10,HR_conv,'LineWidth', lineW);
        % legend([h1,h2],'Global signal', 'X_{HR}')
        % title('BOLD fluctuations due to changes in HR')
        % text(60, 4,  sprintf('r=%3.2f  ',  r_PRF_sc(2)) ,'FontSize',fontTxt,'FontWeight','bold')
        % ylabel('Amplitude (a.u.)')
        % ylim([yl1, yl2])
        % legend('boxoff')
        % 
        % ax3 = subplot(5,3,7:8);
        % h1=plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
        % h2=plot(timeMR,yPred,'LineWidth',lineW);
        % title('Full model')
        % text(60, 4,  sprintf('r=%3.2f  ',  r_PRF_sc(1)) ,'FontSize',fontTxt,'FontWeight','bold')
        % ylabel('Amplitude (a.u.)')
        % legend([h1,h2],'Global signal','X_{FM}')
        % ylim([yl1, yl2])
        % legend('boxoff')
        % 
        % ax4 = subplot(5,3,10:11);
        % h1 = plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
        % h2 = plot(time_10,RF_conv,'LineWidth',lineW);
        % title('BOLD fluctuations due to changes in respiration')
        % text(60, 4,  sprintf('r=%3.2f  ',  r_PRF_sc(3)) ,'FontSize',fontTxt,'FontWeight','bold')
        % legend([h1,h2],'Global signal','X_{RF}'), legend('boxoff')
        % ylabel('Amplitude (a.u.)')
        % ylim([yl1, yl2])
        % 
        % ax7 = subplot(5,3,[12,15]);
        % plot(t_IR,RRF_sc,'LineWidth',4), grid on
        % title('Respiration response function (RRF_{sc}) ')
        % xlim([0 60])
        % xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
        % 
        % ax5 = subplot(5,3,13:14);
        % plot(time_10,RF,'LineWidth',1), hold on
        % title('Respiratory flow (RF)')
        % ylabel('RF (a.u.)')
        % ylim([-0.01 0.1])
        % xlabel('Time (s)')
        % 
        % linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
        % xlim([timeMR(1) timeMR(end)])
        % 
        % ax_list = [ax1,ax2,ax3,ax4,ax5,ax6,ax7];
        % for ax=ax_list
        %     subplot(ax)
        %     ax.XGrid = 'on';
        %     ax.GridAlpha=0.7;
        %     ax.GridLineStyle='--';
        %     ax.FontSize = 17;
        %     ax.FontWeight = 'bold';    
        % end
        
        %%   4: Create matrix of Physiological Regressors for the General linear Model 
        
        xPhys = [HR_conv_MR,RF_conv_MR];  xPhys = detrend(xPhys,'linear');
        
        % figure('Position', [ 316         673        1849         483])
        % plot(timeMR(:), xPhys)
        % xlabel('Time (s)')
        % ylabel('Amplitude (a.u.)')
        % subject = sub_num;
        % title(sprintf('Physiological regressors to be included in the General Linear Model for scan %s (%s) ', subject, task),'Interpreter','none')

        CRF_par = PRF_par(1:4);
        RRF_par = PRF_par(5:8);
        scale_par = PRF_par(9:10);
        save(path_out,'t_IR','CRF_sc','RRF_sc','HR_conv','RF_conv','r_PRF_sc','yPred','xPhys','CRF_par','RRF_par','scale_par')
    end
end







