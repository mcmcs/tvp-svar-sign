function [St_old, DCDhalft_old] = restriction_none(B_old_t, ddelta_t, ggamma_t, Q_t, t)
% checking sign restrictions

% --- Auxiliary function

Dt_old    = diag(exp(ddelta_t/2));
Ct_old    = ggammatoC(ggamma_t);
ssigma_old = (Dt_old*Ct_old*Dt_old);
ssigma_old  = (ssigma_old+ssigma_old')/2;



cond_ssigma_old = cond(ssigma_old);
nv = size(Dt_old,1);
if cond_ssigma_old<10000
    DCDhalft_old = chol(ssigma_old, 'lower');

% 
%     B_tmp  = reshape(B_old_t', [], nv)';
% B1_tmp = B_tmp(:,2:(2+nv-1)); % B1 - coeff on the first lag
% B2_tmp = B_tmp(:,(2+nv):(2+2*nv-1)); % B2 - coeff on the second lag
% 
% matrixF = [B1_tmp B2_tmp;
%      eye(nv),zeros(nv)];
% 
% eigF = max(abs(eig(matrixF)));

    
St_old = 1;%(eigF<1);

else
      St_old=0;
  DCDhalft_old = eye(nv)*1e8;
    
end