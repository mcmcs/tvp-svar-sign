function [shatnew,signew]=kfilter(y,H,F,shat,sig,R,Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL
% y(t)=H*s(t)+e(t)
% s(t)=F*s(t-1)+v(t)
% V(e(t))=R
% V(v(t))=Q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y(:); 
omega=F*sig*F'+Q;
sigma=H*omega*H'+R;

% original
% sigma = (sigma + sigma')/2.0;
k=omega*H'/sigma; 

% % modified
% k=omega*H'*(invChol_mex(sigma)); 

sfor=F*shat; 
ferr=y-H*sfor; 
shatnew=sfor+k*ferr;
signew=omega-k*H*omega;