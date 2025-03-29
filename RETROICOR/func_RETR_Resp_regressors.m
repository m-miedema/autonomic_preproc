function [Regr]=func_RETR_Resp_regressors(resp_f,M,Fs,rsp_phase_interp)
% RETROICOR (Respiratory regressors)

NT=length(resp_f);
Phi=zeros(NT,1);

resp_f = smooth(resp_f, 1*Fs);
resp_der=diff(resp_f); resp_der=[resp_der;resp_der(end)];
resp_der = filloutliers(resp_der,'linear','movmedian',0.5*Fs);

NB=1000;
[Val,edges]=histcounts(resp_f,NB);
sign_prev=1;
for i=1:NT    
    v=resp_f(i);    
    [~,edge]=min(abs(v-edges));
    if (edge-v)>0 
        edge=edge-1; 
    end
    area=sum(Val(1:edge));
    sign_resp=sign(resp_der(i));
    if i > 2 && i < NT-2
        if  abs(sum(abs(resp_der(i-2:i+2)))) < 0.0005
            sign_resp=sign_prev;
        end
        % % check for breath-holds more explicitly
        % if i > round(2*fs) && i < NT-round(2*fs)
        %     if  sum(abs(resp_der(i-round(2*fs):i+round(2*fs)))) < 1
        %         sign_resp=sign(sum(resp_der(i-round(2*fs):i+round(2*fs))));
        %     end
        % end
    end
    sign_prev = sign_resp;
    Phi(i)=pi*area*sign_resp/NT;    
end

% find the places where the sign of the respiratory phase needs to be
% adjusted
% i.e. find each brushed region from QC
search_starts = strfind([0 rsp_phase_interp],[0 1]);
search_ends = strfind([rsp_phase_interp 0],[1 0]);
for s = 1:size(search_starts,2)
    sign_start = sign(Phi(search_starts(s)));
    %Phi(search_starts(s):search_ends(s)) = sign_start*abs(Phi(search_starts(s):search_ends(s)));
    Phi(search_starts(s):search_ends(s)) = linspace(Phi(search_starts(s)),Phi(search_ends(s)),(search_ends(s)-search_starts(s)+1));
end

Regr=zeros(NT,M*2);
for i=1:M
    Regr(:,(i-1)*2+1)=cos(i*Phi);
    Regr(:,i*2)=sin(i*Phi);
end



%%



