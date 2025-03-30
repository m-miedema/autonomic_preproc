% Script to assess intra-class correlation coefficient (ICC)
% for connectivity edges across various task and session conditions.

% Dependent on https://www.mathworks.com/matlabcentral/fileexchange/22099-intraclass-correlation-coefficient-icc
% Figures use https://www.mathworks.com/matlabcentral/fileexchange/74851-daboxplot

close all
clear all

c_methods = {'minimal' 'wPRF'};
c_ICCs_rest = {};
c_ICCs_breathing = {};
c_ICCs_ses02_breathing = {};
c_ICCs_ses03_breathing = {};
c_ICCs_ses02_cold = {};
c_ICCs_intrases03_rest = {};
c_bounds_rest = {};
c_bounds_breathing = {};
c_bounds_ses02 = {};

c_conn_mean_rest = {};
c_conn_mean_breathing = {};
c_conn_mean_ses02 = {};

for c_type = 1:2

    % load in data
    root_dir = "D:\estimation";
    cleaning = c_methods{c_type};
    path_all_conn_output = strcat(root_dir,'\CAN_conn\',cleaning,'_CAN_conn_z.mat');
    load(path_all_conn_output,"all_CAN_conn")

    % load in session 3 intra-scan data as a special case
    path_all_conn_output = strcat(root_dir,'\CAN_conn\',cleaning,'_ses03_rest_CAN_conn_z.mat');
    load(path_all_conn_output,"CAN_conn_a","CAN_conn_b")
    
    % prepare to extract upper triangles
    tri = find(triu(ones(size(all_CAN_conn{1},1,2)),1));
    
    % do comparison of resting state data
    rest02 = all_CAN_conn{1};
    rest03 = all_CAN_conn{4};
    tri02 = [];
    tri03 = [];
    for i = 1:size(rest02,3)
        rest02_i = rest02(:,:,i);
        tri02 = [tri02;rest02_i(tri).'];
        rest03_i = rest03(:,:,i);
        tri03 = [tri03;rest03_i(tri).'];
    end

    c_conn_mean_rest{end+1} = (mean(tri02)+mean(tri03))/2;
    
    ICC_rest = [];
    for j = 1:size(tri02,2)
        M = [tri02(:,j),tri03(:,j)];
        [r, LB, UB, F, df1, df2, p] = ICC(M, 'C-1');%, alpha, r0);
        ICC_rest = [ICC_rest;r,LB,UB];
    end
    c_ICCs_rest{end+1} = ICC_rest(:,1);
    c_bounds_rest{end+1} = ICC_rest(:,2:3);
    
    % do comparison of breathing data
    breathing02 = all_CAN_conn{2};
    breathing03 = all_CAN_conn{5};
    tri02 = [];
    tri03 = [];
    for i = 1:size(breathing02,3)
        breathing02_i = breathing02(:,:,i);
        tri02 = [tri02;breathing02_i(tri).'];
        breathing03_i = breathing03(:,:,i);
        tri03 = [tri03;breathing03_i(tri).'];
    end

    c_conn_mean_breathing{end+1} = (mean(tri02)+mean(tri03))/2;
    
    ICC_breathing = [];
    for j = 1:size(tri02,2)
        M = [tri02(:,j),tri03(:,j)];
        [r, LB, UB, F, df1, df2, p] = ICC(M, 'C-1');%, alpha, r0);
        ICC_breathing = [ICC_breathing;r,LB,UB];
    end

    c_ICCs_breathing{end+1} = ICC_breathing(:,1);
    c_bounds_breathing{end+1} = ICC_breathing(:,2:3);

    % do comparison of resting state and task data
    coldpressor02 = all_CAN_conn{3};
    tri_r = [];
    tri_r2 = [];
    tri_b = [];
    tri_b2 = [];
    tri_c = [];
    for i = 1:size(rest02,3)
        rest02_i = rest02(:,:,i);
        tri_r = [tri_r;rest02_i(tri).'];
        rest03_i = rest03(:,:,i);
        tri_r2 = [tri_r2;rest03_i(tri).'];
        breathing02_i = breathing02(:,:,i);
        tri_b = [tri_b;breathing02_i(tri).'];
        breathing03_i = breathing03(:,:,i);
        tri_b2 = [tri_b2;breathing03_i(tri).'];
        coldpressor02_i = coldpressor02(:,:,i);
        tri_c = [tri_c;coldpressor02_i(tri).'];
    end

    c_conn_mean_ses02{end+1} = (mean(tri_r)+mean(tri_b)+mean(tri_c))/3;
    
    ICC_ses02_breathing = [];
    ICC_ses03_breathing = [];
    for j = 1:size(tri_r,2)
        M = [tri_r(:,j),tri_b(:,j)];
        [r, LB, UB, F, df1, df2, p] = ICC(M, 'C-1');
        ICC_ses02_breathing = [ICC_ses02_breathing;r,LB,UB];

        M = [tri_r2(:,j),tri_b2(:,j)];
        [r, LB, UB, F, df1, df2, p] = ICC(M, 'C-1');
        ICC_ses03_breathing = [ICC_ses03_breathing;r,LB,UB];
    end

    ICC_ses02_cold = [];
    for j = 1:size(tri_r,2)
        M = [tri_r(:,j),tri_c(:,j)];
        [r, LB, UB, F, df1, df2, p] = ICC(M, 'C-1');
        ICC_ses02_cold = [ICC_ses02_cold;r,LB,UB];
    end

    c_ICCs_ses02_breathing{end+1} = ICC_ses02_breathing(:,1);
    c_ICCs_ses03_breathing{end+1} = ICC_ses03_breathing(:,1);
    c_ICCs_ses02_cold{end+1} = ICC_ses02_cold(:,1);
    c_bounds_ses02{end+1} = ICC_ses02_breathing(:,2:3);

    % now do a comparison of the first and second halves of ses03 rest
    tri02 = [];
    tri03 = [];
    for i = 1:size(CAN_conn_a,3)
        CAN_conn_a_i = CAN_conn_a(:,:,i);
        tri02 = [tri02;CAN_conn_a_i(tri).'];
        CAN_conn_b_i = CAN_conn_b(:,:,i);
        tri03 = [tri03;CAN_conn_b_i(tri).'];
    end
    % calculate ICC
    ICC_intra_rest = [];
    for j = 1:size(tri_r,2)
        M = [tri02(:,j),tri03(:,j)];
        [r, LB, UB, F, df1, df2, p] = ICC(M, 'C-1');
        ICC_intra_rest = [ICC_intra_rest;r,LB,UB];
    end

    c_ICCs_intrases03_rest{end+1} = ICC_intra_rest(:,1);
end

%% make box plots
% and make a scatter plot correlating ICC and mean connectivity value
path_boxfig = strcat(root_dir,'\CAN_conn\ICC_boxplots_z.png');
boxfig = figure('units','normalized','outerposition',[0 0 1 1]);
t1 = tiledlayout(1,4,'TileSpacing','compact');
nexttile
h = daboxplot(c_ICCs_rest,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', c_methods); 
title('ICCs for resting state connectivity edges')
ylabel('ICC value')
nexttile
h = daboxplot(c_ICCs_breathing,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', c_methods); 
title('ICCs for breathing connectivity edges')
ylabel('ICC value')
nexttile

c_ICCs_both_breathing = {[c_ICCs_ses02_breathing{1},c_ICCs_ses03_breathing{1}],[c_ICCs_ses02_breathing{2},c_ICCs_ses03_breathing{2}]};
h = daboxplot(c_ICCs_both_breathing,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', c_methods, 'legend',{'session 1', 'session 2'}); 
title('ICCs for inter-condition connectivity edges (rest vs. breathing, sesssions 1 & 2)')
ylabel('ICC value')
nexttile
h = daboxplot(c_ICCs_ses02_cold,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels', c_methods); 
title('ICCs for inter-condition connectivity edges (rest vs. cold pressor, session 1)')
ylabel('ICC value')
saveas(boxfig,path_boxfig)

figure
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
c_ICCs_minimal = [c_ICCs_rest{1}, c_ICCs_breathing{1}, ...
    c_ICCs_ses02_breathing{1},c_ICCs_ses03_breathing{1},c_ICCs_ses02_cold{1},...
    c_ICCs_intrases03_rest{1}];
c_ICCs_wPRF = [c_ICCs_rest{2}, c_ICCs_breathing{2}, ...
    c_ICCs_ses02_breathing{2},c_ICCs_ses03_breathing{2},c_ICCs_ses02_cold{2},...
    c_ICCs_intrases03_rest{2}];
c_ICCs_all = {c_ICCs_minimal,c_ICCs_wPRF};
comp_names = {'rest (session 1 vs. 2)','breathing (session 1 vs. 2)',...
    'rest vs. breathing (session 1)','rest vs. breathing (session 2)',...
    'rest vs. cold pressor (session 1)','rest (session 2 split halves)'};
yline(0.4,':b')
hold on
yline(0.6,':b')
h = daboxplot(c_ICCs_all,'scatter',2,'whiskers',0,'boxalpha',0.7,...
    'xtlabels',comp_names, 'legend', c_methods,'color',c(1:2:3,:)); 
ylim([-0.6 0.9])
xL = xlim;
yL = ylim;
xlim([xL(1) xL(2)+0.5]);
text(xL(2)+0.5,0.4,'acceptable','Color',[0.2 0.4 0.8],'FontSize',8,...
    'HorizontalAlignment','right','VerticalAlignment','cap')
text(xL(2)+0.5,0.6,'good','Color',[0.2 0.4 0.8],'FontSize',8,...
    'HorizontalAlignment','right','VerticalAlignment','cap')
title('ICCs between sessions and conditions')

path_scatterfig = strcat(root_dir,'\CAN_conn\ICCxConn_z.png');
lwidth=1.4;
scatterfig = figure('units','normalized','outerposition',[0 0 1 1]);
t1 = tiledlayout(1,3,'TileSpacing','compact');
nexttile
plot(c_conn_mean_rest{1},c_ICCs_rest{1},'o','LineWidth',lwidth)
hold on
plot(c_conn_mean_rest{2},c_ICCs_rest{2},'o','LineWidth',lwidth)
legend({'Rest, minimal' 'Rest, with PRF'},'Location','southeast')
xlabel('Mean connectivity across subjects')
ylabel('ICC score')
nexttile
plot(c_conn_mean_breathing{1},c_ICCs_breathing{1},'^','LineWidth',lwidth)
hold on
plot(c_conn_mean_breathing{2},c_ICCs_breathing{2},'^','LineWidth',lwidth)
legend({'Breathing, minimal' 'Breathing, with PRF'},'Location','southeast')
xlabel('Mean connectivity across subjects')
ylabel('ICC score')
nexttile
plot(c_conn_mean_rest{1},c_ICCs_rest{2}-c_ICCs_rest{1},'o','Color',[0.4660 0.6740 0.1880],'LineWidth',lwidth)
hold on
plot(c_conn_mean_breathing{1},c_ICCs_breathing{2}-c_ICCs_breathing{1},'^','Color',[0.3660 0.8740 0.6880],'LineWidth',lwidth)
ylabel('Change in ICC score')
xlabel('Mean connectivity across subjects')
legend({'Rest' 'Breathing'},'Location','southeast')
saveas(scatterfig,path_scatterfig)

%% plot connectivity maps with ICC scores overlaid
breathingPRFfig = figure;
path_breathingPRFfig = strcat(root_dir,'\CAN_conn\conn_ICC_breathing_wPRF_z.png');
[colormap2]=cbrewer2('seq','OrRd',100);
imshow(mean(all_CAN_conn{2},3),[0 0.75],'InitialMagnification','fit')%[-0.2 0.2],
hold on
% plot ICC values on top
% text(100,100,[ '[' num2str(x), num2str(y) ']'])
ICC_str = strings(size(all_CAN_conn{2},1,2));
ICC_str(tri) = num2str(round(c_ICCs_breathing{2},2));
[X,Y] = meshgrid(1:size(ICC_str,1), 1:size(ICC_str,1));
text(X(:)-0.45,Y(:),ICC_str(:),'FontSize',5)
xlabel('ROIs in the CAN')
ax = gca;
ax.FontSize = 7; 
colormap(colormap2)
cb = colorbar;
cb.Label.String = 'Mean correlation for session 2 breathing';
cb.Label.FontAngle = 'italic';
title('Intraclass Correlation Coefficients for breathing task, with PRF')
saveas(breathingPRFfig,path_breathingPRFfig)

restPRFfig = figure;
path_restPRFfig = strcat(root_dir,'\CAN_conn\conn_ICC_rest_wPRF_z.png');
[colormap2]=cbrewer2('seq','OrRd',100);
imshow(mean(all_CAN_conn{1},3),[0 0.75],'InitialMagnification','fit')%[-0.2 0.2],
hold on
% plot ICC values on top
% text(100,100,[ '[' num2str(x), num2str(y) ']'])
ICC_str = strings(size(all_CAN_conn{1},1,2));
ICC_str(tri) = num2str(round(c_ICCs_rest{2},2));
[X,Y] = meshgrid(1:size(ICC_str,1), 1:size(ICC_str,1));
%text(19,19,'test')
text(X(:)-0.45,Y(:),ICC_str(:),'FontSize',5)
xlabel('ROIs in the CAN')
ax = gca;
ax.FontSize = 7; 
colormap(colormap2)
cb = colorbar;
cb.Label.String = 'Mean correlation for session 2 rest';
cb.Label.FontAngle = 'italic';
title('Intraclass Correlation Coefficients for rest, with PRF')
saveas(restPRFfig,path_restPRFfig)


% NOW RE-LOAD FOR MINIMAL AGAIN
cleaning = 'minimal';
path_all_conn_output = strcat(root_dir,'\CAN_conn\',cleaning,'_CAN_conn_z.mat');
load(path_all_conn_output,"all_CAN_conn")

breathingminfig = figure;
path_breathingminfig = strcat(root_dir,'\CAN_conn\conn_ICC_breathing_minimal_z.png');
[colormap2]=cbrewer2('seq','OrRd',100);
imshow(mean(all_CAN_conn{2},3),[0 0.75],'InitialMagnification','fit')%[-0.2 0.2],
hold on
% plot ICC values on top
% text(100,100,[ '[' num2str(x), num2str(y) ']'])
ICC_str = strings(size(all_CAN_conn{2},1,2));
ICC_str(tri) = num2str(round(c_ICCs_breathing{1},2));
[X,Y] = meshgrid(1:size(ICC_str,1), 1:size(ICC_str,1));
%text(19,19,'test')
text(X(:)-0.45,Y(:),ICC_str(:),'FontSize',5)
xlabel('ROIs in the CAN')
ax = gca;
ax.FontSize = 7; 
colormap(colormap2)
cb = colorbar;
cb.Label.String = 'Mean correlation for session 2 breathing';
cb.Label.FontAngle = 'italic';
title('Intraclass Correlation Coefficients for breathing task, minimal')
saveas(breathingminfig,path_breathingminfig)

restminfig = figure;
path_restminfig = strcat(root_dir,'\CAN_conn\conn_ICC_rest_minimal_z.png');
[colormap2]=cbrewer2('seq','OrRd',100);
imshow(mean(all_CAN_conn{1},3),[0 0.75],'InitialMagnification','fit')%[-0.2 0.2],
hold on
% plot ICC values on top
% text(100,100,[ '[' num2str(x), num2str(y) ']'])
ICC_str = strings(size(all_CAN_conn{1},1,2));
ICC_str(tri) = num2str(round(c_ICCs_rest{1},2));
[X,Y] = meshgrid(1:size(ICC_str,1), 1:size(ICC_str,1));
%text(19,19,'test')
text(X(:)-0.45,Y(:),ICC_str(:),'FontSize',5)
xlabel('ROIs in the CAN')
ax = gca;
ax.FontSize = 7; 
colormap(colormap2)
cb = colorbar;
cb.Label.String = 'Mean correlation for session 2 rest';
cb.Label.FontAngle = 'italic';
title('Intraclass Correlation Coefficients for rest, minimal')
saveas(restminfig,path_restminfig)

% output some stats
[h,p] = ttest(c_ICCs_rest{1},c_ICCs_rest{2})
[h,p] = ttest(c_ICCs_breathing{1},c_ICCs_breathing{2})
[h,p] = ttest(c_ICCs_ses02_breathing{1},c_ICCs_ses02_breathing{2})

% re-output some stats
for i = 1:5
    disp(comp_names(i));
    [h,p] = ttest(c_ICCs_minimal(:,i),c_ICCs_wPRF(:,i))
end

% apply an ANOVA to test whether within-session significantly different
% from inter-task
aov = anova(c_ICCs_wPRF);
aov = anova(c_ICCs_minimal);

% factor 1: pipeline
f_pips = [ones(size(c_ICCs_minimal));2*ones(size(c_ICCs_wPRF))];
% factor 2: specific task types being examined
f_tasks = [ones(171,1),2*ones(171,1),3*ones(171,2),4*ones(171,1);...
    ones(171,1),2*ones(171,1),3*ones(171,2),4*ones(171,1)];
% factor 3: edge (which ROI to which ROI)
ICC_ind = 1:171;
ICC_ind = ICC_ind.';
f_ICCs = repmat(ICC_ind,[2 5]);

y = [c_ICCs_minimal;c_ICCs_wPRF];
[p,tbl,stats] = anovan(y(:),{f_pips(:),f_tasks(:),f_ICCs(:)},'model','interaction','varnames',{'pipeline','condition','edge'});

writecell(tbl,'anovaresults.xlsx')