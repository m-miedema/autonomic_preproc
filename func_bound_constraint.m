function [c, ceq] = bound_constraint(x)
%SIMPLE_CONSTRAINT Nonlinear inequality constraints for PRF estimation.

tol = 0.01; % tolerance for limit values of PRF
t_edge = 60; % time when the PRF should approach this limit

t_win= 0:0.1:t_edge;

% calculate the PRF values at the edge
t1c=x(1) ; d1c=x(2);
t2c=x(3);  d2c=x(4);
t1r=x(5); d1r = x(6);
t2r=x(7); d2r=x(8);
R_CRF = x(9); R_RRF = x(10);

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
C_edge = CRF_sc(end);
R_edge = RRF_sc(end);

% a1= sqrt(t1c)/d1c; a2= sqrt(t1c)*d1c ;
% IR_cardF = (t_edge^a1)*exp(-t_edge/a2);  
% a1= sqrt(t2c)/d2c; a2= sqrt(t2c)*d2c ;
% IR_cardS = (t_edge^a1)*exp(-t_edge/a2);     
% 
% a1= sqrt(t1r)/d1r; a2= sqrt(t1r)*d1r ;
% IR_respF = (t_edge^a1)*exp(-t_edge/a2);  
% a1= sqrt(t2r)/d2r; a2= sqrt(t2r)*d2r ;
% IR_respS = (t_edge^a1)*exp(-t_edge/a2);  

%CRF_edge = IR_cardF + R_CRF*IR_cardS;     
%RRF_edge = IR_respF + R_RRF*IR_respS;  
% note that the scaling factors of the PRFs are being neglected here -
% might need to calculate this for overall window for this to work

% dbstop in bound_constraint at 29 if ~isreal(IR_cardF)
% dbstop in bound_constraint at 29 if ~isreal(IR_cardS)
% dbstop in bound_constraint at 29 if ~isreal(IR_respF)
% dbstop in bound_constraint at 29 if ~isreal(IR_respS)
% 
% if isnan(IR_cardF)
%     IR_cardF = 1;
% end
% if isnan(IR_cardS)
%     IR_cardS = 1;
% end
% if isnan(IR_respF)
%     IR_respF = 1;
% end
% if isnan(IR_respS)
%     IR_respS = 1;
% end
if isnan(C_edge)
    C_edge = 1;
end
if isnan(R_edge)
    R_edge = 1;
end

c = [abs(C_edge)-tol; abs(R_edge)-tol];
% c(1) = abs(IR_cardS)-tol;
% c(2) = abs(IR_cardF)-tol;
% c(3) = abs(IR_respS)-tol;
% c(4) = abs(IR_respF)-tol;

% no nonlinear equality constraints:
ceq = [];