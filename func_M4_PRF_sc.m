function [obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred,yPred_card,yPred_resp,HR_conv_MR, RF_conv_MR,B]=func_M4_PRF_sc(P,Ts_10,HR,RF,ind_BOLD_10,GS,only_OF)
% I changed this to add the yPred variants without updating
% estimate_scan_PRF
t1c=P(1) ; d1c=P(2);
t2c=P(3);  d2c=P(4);
t1r=P(5); d1r = P(6);
t2r=P(7); d2r=P(8);
R_CRF = P(9); R_RRF = P(10);

NV = length(GS);
t_win= 0 :Ts_10:60;

a1= sqrt(t1c)/d1c; a2= sqrt(t1c)*d1c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardF=IR/max(IR);  
a1= sqrt(t2c)/d2c; a2= sqrt(t2c)*d2c ;
IR = t_win.^a1.*exp(-t_win/a2); IR_cardS=IR/max(IR);     

a1= sqrt(t1r)/d1r; a2= sqrt(t1r)*d1r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respF=IR/max(IR);   
a1= sqrt(t2r)/d2r; a2= sqrt(t2r)*d2r ;
IR = t_win.^a1.*exp(-t_win/a2); IR_respS=IR/max(IR);  

CRF_sc = IR_cardF + R_CRF*IR_cardS;     
RRF_sc = IR_respF + R_RRF*IR_respS;  
CRF_sc = CRF_sc/max(abs(CRF_sc));
RRF_sc = RRF_sc/max(abs(RRF_sc));

r_PRF_sc=zeros(3,1);

HR = zscore(HR); RF = zscore(RF); GS = zscore(GS);
% HR_Fconv=conv(HR,IR_cardF); HR_Fconv_MR=HR_Fconv(ind_BOLD_10); % the convolution is causing the errors here
% HR_Sconv=conv(HR,IR_cardS); HR_Sconv_MR=HR_Sconv(ind_BOLD_10);
% RF_Fconv=conv(RF,IR_respF); RF_Fconv_MR=RF_Fconv(ind_BOLD_10);
% RF_Sconv=conv(RF,IR_respS); RF_Sconv_MR=RF_Sconv(ind_BOLD_10);
% regr = [HR_Fconv_MR,HR_Sconv_MR,RF_Fconv_MR,RF_Sconv_MR];  regr = detrend(regr,'linear');
% regr = [ones(NV,1),regr];

HR_conv = conv(HR,CRF_sc); HR_conv_MR = HR_conv(ind_BOLD_10);
RF_conv = conv(RF,RRF_sc); RF_conv_MR = RF_conv(ind_BOLD_10); 
regr = [HR_conv_MR,RF_conv_MR];     regr = detrend(regr,'linear');
regr = [ones(NV,1), regr];

B = regr\GS;     yPred = regr*B;
obj_function = 1 - corr(yPred,GS) ;

if only_OF == 0
    % CRF_sc = B(2) * IR_cardF + B(3) * IR_cardS;     CRF_sc = CRF_sc/max(abs(CRF_sc));
    % RRF_sc = B(4)*IR_respF + B(5)*IR_respS; RRF_sc = RRF_sc/max(abs(RRF_sc));
    % 
    % HR_conv = B(2)*HR_Fconv + B(3)*HR_Sconv;    HR_conv = HR_conv(1:length(HR));   HR_conv_MR = HR_conv(ind_BOLD_10);
    % RF_conv = B(4)*RF_Fconv + B(5)*RF_Sconv;     RF_conv = RF_conv(1:length(RF));     RF_conv_MR = RF_conv(ind_BOLD_10);
    % 
    r_PRF_sc(1) = corr(yPred,GS);
    yPred_card = regr(:,2)*B(2);  r_PRF_sc(2) = corr(yPred_card,GS);
    yPred_resp = regr(:,3)*B(3);  r_PRF_sc(3) = corr(yPred_resp,GS);

end



%%   ----------------------------------------------------