
clear, close all, clc
%%  1:  Load physiological variables (heart rate and respiration) and global signal (GS) from MAT-File

%% set up loop across subjects and choose which scan to estimate the group PRF for
subjects = {'0059','0704','0134','1092','1687','1693','1750','3046','4901','5262','5571','5931',...
    '6293','7002','7712','7994','8116','8223','8420','8803','8907','9922'};
nScans = length(subjects);

ses_num="02";
task="breathing";
NV = 397;

% define path to root directory
root_dir = "D:\estimation";
physio_out_dir = "D:\estimation\physio_preproc\";
path_out = strcat(root_dir,'\PRF_group\',task,'_',ses_num,'_PRF.mat');

%% check if any subjects have less physio data than others
n_physio = zeros(nScans,1);
t_zeros = zeros(nScans,1);
for sub_i = 1:nScans
    sub_num=subjects{sub_i};

    ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
    bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/",ext_root);       
    path_physio = strcat(physio_out_dir,bids_root,'_physio_and_triggers.mat');

    % get physio signals
    load(path_physio,"time_10")
    n_physio(sub_i) = length(time_10);
    [~,t_zeros(sub_i)] = min(abs(time_10));
end
n_opt = median(n_physio);
zer_opt = median(t_zeros);
GS_all = zeros(NV,nScans);
HR_all = zeros(n_opt,nScans);
RF_all = zeros(n_opt,nScans);

%% load in the data
for sub_i = 1:length(subjects)
    sub_num=subjects{sub_i};

    ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
    bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/",ext_root);
       
    path_physio = strcat(physio_out_dir,bids_root,'_physio_and_triggers.mat');
    path_GS = strcat(root_dir,'\extracted\',ext_root,'_GS.mat');        
    
    % get global signal        
    load(path_GS,"WB_clean")
    GS = WB_clean; % test this

    % get physio signals
    load(path_physio,"HR_10","resp_10","RF","trig_ind_10","time_10")

    %% line everything up
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
    if size(RF,2) > 1
        RF = RF.';
    end

    n_i = length(HR_10);
    % check if data was too short here
    if n_i < n_opt
        % need to pad the lengths of things accordingly
        extend_pre = t_zeros(sub_i) - zer_opt;
        extend_post = (n_opt - n_i) - extend_pre;
        HR_10 = [HR_10(1)*ones(extend_pre,1);HR_10;HR_10(end)*ones(extend_post,1)]; 
        RF = [RF(1)*ones(extend_pre,1);RF;RF(end)*ones(extend_post,1)]; 
    elseif n_i > n_opt
        warning('Warning: issue with data length, check closely before continuing.')
        disp(sub_num)
        extend_pre = t_zeros(sub_i) - zer_opt;
        extend_post = (n_opt - n_i) - extend_pre;
        HR_10 = HR_10(extend_pre+1:end+extend_post);
        RF = RF(extend_pre+1:end+extend_post);
    end

    GS_all(:,sub_i) = GS;
    HR_all(:,sub_i) = HR_10;
    RF_all(:,sub_i) = RF;
end

% make sure we're using a subject of appropriate length here

ind_BOLD_10 = trig_ind_10;

Ts_10 = 0.1 ;                                                       % Sampling period in seconds
% time_10 = 0:Ts_10:(NP-1)*Ts_10;
% timeMR = time_10(ind_BOLD_10);


%% 2: Estimate PRF parameters

ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',100,'Display','iter','UseParallel',1);   % Display: iter
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'MaxIterations',200,'MaxFunctionEvaluations',3000 ,'OptimalityTolerance',1e-8 );

PRF_par = [  3.1    2.5   5.6    0.9    1.9   2.9   12.5    0.5  ];  
lb(1:8)=0;  ub(1:8)=20; ub(2:2:8)=3;   lb(9:10) = -inf; ub(9:10)=inf;

h = @(P) func_M3_PRF_pop_sc(P,Ts_10,HR_all,RF_all,ind_BOLD_10,GS_all);

tic
% Uncomment the following line if you want to use  Genetic Algorithm
% (GA). GA may yield better fit with the cost of longer computational time.
% PRF_par = ga(h,length(lb),[],[],[],[],lb,ub,[],[],ga_opts);
PRF_par = fmincon(h,PRF_par,[],[],[],[],lb,ub,[],options);
fprintf('Minutes elapsed: %3.1f  \n',  toc/60)


 [obj_function,CRF_1, CRF_2, RRF_1, RRF_2 ,r_PRF_pop ] = h(PRF_par);

fprintf(' ----------------------------------------------- \n')
fprintf('Correlation b/w GS and PRF output \n')
fprintf('CRF & RRF (HR & RF): r=%3.2f  \n',mean(r_PRF_pop))

%%  3: Plot PRF curves  

t_IR = 0:Ts_10:(length(CRF_1)-1)*Ts_10;

figure
subplot(2,1,1)
plot(t_IR,CRF_1,'LineWidth',2), grid on, hold on
plot(t_IR,CRF_2,'LineWidth',2)
title('Cardiac Response Function (CRF_{pop}) ')
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
xlim([0 60])

subplot(2,1,2)
plot(t_IR,RRF_1,'LineWidth',2), grid on, hold on
plot(t_IR,RRF_2,'LineWidth',2)
title('Respiration response function (RRF_{pop}) ')
xlim([0 60])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

%%  4: Apply PRF model (timeseries and PRF curves)  

%  Set the following parameters !!
sub_i = 20;
sc = sub_i;     % choose a scan (sc) from 1-164

% -----------------------------------------

subject = subjects{sub_i};
fprintf('Scan: %d,    Subject: %s,    Task: %s    \n', sc, subject, task)

GS = zscore(GS_all(:,sc)); HR = zscore(HR_all(:,sc)); RF = zscore(RF_all(:,sc));
NV = length(GS);

HR_conv_1 = conv(HR,CRF_1); HR_conv_1_MR = HR_conv_1(ind_BOLD_10);
HR_conv_2 = conv(HR,CRF_2); HR_conv_2_MR = HR_conv_2(ind_BOLD_10);
RF_conv_1 = conv(RF,RRF_1); RF_conv_1_MR = RF_conv_1(ind_BOLD_10);
RF_conv_2 = conv(RF,RRF_2); RF_conv_2_MR = RF_conv_2(ind_BOLD_10);

xPhys = [HR_conv_1_MR, HR_conv_2_MR, RF_conv_1_MR, RF_conv_2_MR];     xPhys = detrend(xPhys,'linear');
regr = [ones(NV,1), xPhys];

B = regr\GS;     yPred = regr*B;
r_PRF(1) = corr(GS, yPred);

yPred_card = regr(:,2:3)*B(2:3);  r_PRF(2) = corr(yPred_card,GS);
yPred_resp = regr(:,4:5)*B(4:5);  r_PRF(3) = corr(yPred_resp,GS);

CRF = B(2) * CRF_1 + B(3) * CRF_2;     CRF = CRF/max(abs(CRF));
RRF = B(4)*RRF_1 + B(5)*RRF_2; RRF = RRF/max(abs(RRF));

HR_conv = B(2)*HR_conv_1 + B(3)*HR_conv_2;    HR_conv = HR_conv(1:length(HR));   HR_conv_MR = HR_conv(ind_BOLD_10);
RF_conv = B(4)*RF_conv_1 + B(5)*RF_conv_2;     RF_conv = RF_conv(1:length(RF));     RF_conv_MR = RF_conv(ind_BOLD_10);


fprintf(' ----------------------------------------------- \n')
fprintf('Correlation b/w GS and PRF output \n')
fprintf('CRF (HR): %3.2f  \n',r_PRF(2))
fprintf('RRF (RF): %3.2f  \n',r_PRF(3))
fprintf('CRF & RRF (HR & RF): %3.2f  \n',r_PRF(1))


%%  5: Plot output of PRF model (timeseries and PRF curves)  

%  Set the following parameters !!

smoothPar = 5;
fontTitle = 20;
fontLabels = 8;
fontTxt = 12;
lineW = 3;
yl1 = -5.3; yl2 = 5.5;

% -----------------------------------------

HR = HR_all(:,sc);
RF = RF_all(:,sc);

figure

ax1 = subplot(5,3,1:2);
plot(time_10,HR)
ylabel('HR (bpm)')
title(sprintf('Heart rate (HR; %2.0f±%1.0f bpm )',mean(HR),std(HR)))

ax6 = subplot(5,3,[3,6]);
plot(t_IR,CRF,'LineWidth',4), grid on
title('Cardiac Response Function (CRF_{sc}) ')
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')
xlim([0 60])

ax2 = subplot(5,3,4:5);
h1=plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
h2=plot(time_10,HR_conv,'LineWidth', lineW);
legend([h1,h2],'Global signal', 'X_{HR}')
title('BOLD fluctuations due to changes in HR')
text(60, 4,  sprintf('r=%3.2f  ',  r_PRF(2)) ,'FontSize',fontTxt,'FontWeight','bold')
ylabel('Amplitude (a.u.)')
ylim([yl1, yl2])
legend('boxoff')


ax3 = subplot(5,3,7:8);
h1=plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
h2=plot(timeMR,yPred,'LineWidth',lineW);
title('Full model')
text(60, 4,  sprintf('r=%3.2f  ',  r_PRF(1)) ,'FontSize',fontTxt,'FontWeight','bold')
ylabel('Amplitude (a.u.)')
legend([h1,h2],'Global signal','X_{FM}')
ylim([yl1, yl2])
legend('boxoff')

ax4 = subplot(5,3,10:11);
h1 = plot(timeMR,smooth(GS,smoothPar),'LineWidth',lineW); hold on
h2 = plot(time_10,RF_conv,'LineWidth',lineW);
title('BOLD fluctuations due to changes in respiration')
text(60, 4,  sprintf('r=%3.2f  ',  r_PRF(3)) ,'FontSize',fontTxt,'FontWeight','bold')
legend([h1,h2],'Global signal','X_{RF}'), legend('boxoff')
ylabel('Amplitude (a.u.)')
ylim([yl1, yl2])


ax7 = subplot(5,3,[12,15]);
plot(t_IR,RRF,'LineWidth',4), grid on
title('Respiration response function (RRF_{sc}) ')
xlim([0 60])
xlabel('Time (s)'), ylabel('Amplitude (a.u.)')

ax5 = subplot(5,3,13:14);
plot(time_10,RF,'LineWidth',1), hold on
title('Respiratory flow (RF)')
ylabel('RF (a.u.)')
xlabel('Time (s)')

linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
xlim([timeMR(1) timeMR(end)])

ax_list = [ax1,ax2,ax3,ax4,ax5,ax6,ax7];
for ax=ax_list
    subplot(ax)
    ax.XGrid = 'on';
    ax.GridAlpha=0.7;
    ax.GridLineStyle='--';
    ax.FontSize = 13;
    ax.FontWeight = 'bold';    
end

%%   6: Create matrix of Physiological Regressors for the General linear Model 

% figure('Position', [ 316         673        1849         483])
% 
% plot(timeMR, xPhys)
% xlabel('Time (s)')
% ylabel('Amplitude (a.u.)')

subject = subjects{sub_i};
title(sprintf('Physiological regressors to be included in the General Linear Model for scan %s (%s) ', subject, task),'Interpreter','none')











