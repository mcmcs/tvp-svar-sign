function [St_old, DCDhalft_old] = restriction_simple3(B_t, ddelta_t, ggamma_t, Q_t, t)
% simple sign restriction function
% this is for DGP to do Geweke test

% NOTE: B_t is bbet_t in our paper

% --- first restriction
Dt_old    = diag(exp(ddelta_t/2));
Ct_old    = ggammatoC(ggamma_t);
DCDhalft_old = chol(Dt_old*Ct_old*Dt_old, 'lower');
% Lt_old = DCDhalft_old * Q_t;
% St_old_1 = (Lt_old(1,1) > 0) & (Lt_old(2,1) < 0);

St_old_1 = true;

% --- second restriction 

St_old_2 = true;

% % Restriction used in the paper
% % ---
% eqn_ffr = 1; % 1 because we are identifying the first shock
% ind_ffr = 3; % 3 because ffr is in the third location in the variable ordering
% ind_y  = 1; % 2 because u is in the first location in the variable ordering
% ind_p  = 2; % 1 because infl is in the first location in the variable ordering
% ind_m  = 4;
% ind_cs = 5;
% hrz_compute = 60; % horizon to compute
% hrz_impose  = 60; % horizon to impose
% % ---
% 
% nv = size(Dt_old,1);
% B_tmp  = reshape(B_t', [], nv)';
% B1_tmp = B_tmp(:,2:(2+nv-1)); % B1 - coeff on the first lag
% B2_tmp = B_tmp(:,(2+nv):(2+2*nv-1)); % B2 - coeff on the second lag
% 
% Ap_tmp = inv(Lt_old); %A(t)'
% F1_tmp = (Ap_tmp*B1_tmp); % F1 - coeff on the first lag
% F2_tmp = (Ap_tmp*B2_tmp); % F2 - coeff on the second lag
% res_to_p = impulse_to_pmt(Ap_tmp, F1_tmp, F2_tmp, hrz_compute, eqn_ffr, ind_ffr, ind_p); % response of ffr to permanent chg in pi
% 
% St_old_2 = (res_to_p(hrz_impose)<4) * (res_to_p(hrz_impose)> 1); % taylor principle

% % For future (stationary restriction)
% nv = size(Dt_old,1);
% B_tmp  = reshape(B_t', [], nv)';
% B1_tmp = B_tmp(:,2:(2+nv-1)); % B1 - coeff on the first lag
% B2_tmp = B_tmp(:,(2+nv):(2+2*nv-1)); % B2 - coeff on the second lag
% B_comp = [B1_tmp, B2_tmp; ...
%  eye(nv), zeros(nv)]; %companion form
% St_old_2 = max ( abs( eig(B_comp) ) ) < 1; %stationary restriction

% --- final

St_old = St_old_1 & St_old_2;
