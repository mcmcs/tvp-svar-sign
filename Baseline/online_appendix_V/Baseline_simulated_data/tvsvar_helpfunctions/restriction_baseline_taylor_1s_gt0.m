function [St_old, DCDhalft_old] = restriction_baseline_taylor_1s_gt0(B_old_t, ddelta_t, ggamma_t, Q_t, t)
% checking sign restrictions

% _lg: larger bound
% _1s: one-sided >1

% ===============
% Information about the sign restrictions
% systematic component identification

% --- Because we have time dependent restrictions
ttime=1969.75:0.25:2023.25;
% ttime_zlb = [1979.75:0.25:1982.75,2009:0.25:2015.5,2020.25:0.25:2021.75];
% ttime_ffr = setdiff(ttime,ttime_zlb);
ttime_ffr = [1969.75:0.25:1979.50, ...
    1983.00:0.25:2008.75, ...
    2015.75:0.25:2020.00, ...
    2022.00:0.25:2023.25];

ttime_zlb0 = [1979.75:0.25:1982.75];
ttime_zlb1 = [2009:0.25:2015.5];
ttime_zlb2 = [2020.25:0.25:2021.75];

% --- Extra info for restrictions
eqn_ffr = 1; % 1 because we are identifying the first shock
ind_ffr = 3; % 3 because ffr is in the third location in the variable ordering
ind_y  = 1; % 2 because u is in the first location in the variable ordering
ind_p  = 2; % 1 because infl is in the first location in the variable ordering
ind_m  = 4;
ind_cs = 5;
hrz_compute = 60; % horizon to compute
hrz_impose  = 60; % horizon to impose
% hrz_compute = Inf; % horizon to compute
% hrz_impose  = 1; % horizon to impose
% ===============

% --- Auxiliary function

Dt_old    = diag(exp(ddelta_t/2));
Ct_old    = ggammatoC(ggamma_t);
ssigma_old = (Dt_old*Ct_old*Dt_old);
ssigma_old  = (ssigma_old+ssigma_old')/2;

cond_ssigma_old = cond(ssigma_old);
nv = size(Dt_old,1);
if cond_ssigma_old<10000

    DCDhalft_old = chol(ssigma_old, 'lower');

    % B_tmp  = reshape(B_old_t', [], nv)';
    % B1_tmp = B_tmp(:,2:(2+nv-1)); % B1 - coeff on the first lag
    % B2_tmp = B_tmp(:,(2+nv):(2+2*nv-1)); % B2 - coeff on the second lag
    %
    % matrixF = [B1_tmp B2_tmp;
    %      eye(nv),zeros(nv)];
    %
    % eigF = max(abs(eig(matrixF)));

    % --- Check sign restrictions
    L0 = DCDhalft_old * Q_t;
    At = inv(L0)';
    ppsi_y = -At(1,1)/At(3,1);
    ppsi_p = -At(2,1)/At(3,1);
    ppsi_m = -At(4,1)/At(3,1);
    ppsi_cs = -At(5,1)/At(3,1);

    % --- effect of permanent changes
    nv = size(Dt_old,1);
    B_tmp  = reshape(B_old_t', [], nv)';
    B1_tmp = B_tmp(:,2:(2+nv-1)); % B1 - coeff on the first lag
    B2_tmp = B_tmp(:,(2+nv):(2+2*nv-1)); % B2 - coeff on the second lag
    Ap_tmp = At';
    F1_tmp = (Ap_tmp*B1_tmp); % F1 - coeff on the first lag
    F2_tmp = (Ap_tmp*B2_tmp); % F2 - coeff on the second lag
    % res_to_y  = impulse_to_pmt(Ap_tmp, F1_tmp, F2_tmp, hrz_compute, eqn_ffr, ind_ffr, ind_y); % response of ffr to permanent chg in u
    res_to_p = impulse_to_pmt(Ap_tmp, F1_tmp, F2_tmp, hrz_compute, eqn_ffr, ind_ffr, ind_p); % response of ffr to permanent chg in pi
    % res_to_m  = impulse_to_pmt(Ap_tmp, F1_tmp, F2_tmp, hrz_compute, eqn_ffr, ind_ffr, ind_m); % response of ffr to permanent chg in pi
    % res_to_cs = impulse_to_pmt(Ap_tmp, F1_tmp, F2_tmp, hrz_compute, eqn_ffr, ind_ffr, ind_cs); % response of ffr to permanent chg in pi

    % --- Actual checks
    % sign_indicator_ffr = (At(1,1)<0)*(At(2,1)<0)*(At(3,1)>0)*(At(4,1)<0)*(At(5,1)>0)*(L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0)*(ppsi_y<4)*(ppsi_p<4)*(ppsi_m<4)*(ppsi_cs>-4);
    sign_indicator_zlb0  = (L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0);
    sign_indicator_zlb1  = (L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0);
    sign_indicator_zlb2  = (L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0);

    % sign_indicator_ffr1 = (At(1,1)<0)*(At(2,1)<0)*(At(3,1)>0)*(At(4,1)<0)*(At(5,1)>0)*(L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0)*(L0(5,1)>0)*(ppsi_y<4)*(ppsi_p<4)*(ppsi_m<4)*(ppsi_cs>-4);
    % sign_indicator_ffr1 = (At(1,1)<0)*(At(3,1)>0)*(At(4,1)<0)*(At(5,1)>0)*(L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0)*(L0(5,1)>0)*(ppsi_y<4)*(ppsi_m<4)*(ppsi_cs>-4); % replace original with taylor

    sign_indicator_ffr1 = (At(1,1)<0)*(At(2,1)<0)*(At(3,1)>0)*(At(4,1)<0)*(At(5,1)>0)*(L0(3,1)>0)*(L0(4,1)<0)*(L0(2,1)<0)*(ppsi_y<4)*(ppsi_p<4)*(ppsi_m<4)*(ppsi_cs>-4);
    sign_indicator_ffr2 = (res_to_p(hrz_impose)>0); % taylor principle
    sign_indicator_ffr = sign_indicator_ffr1 & sign_indicator_ffr2; %combine both

    St_old = ((sum(ttime(t)==ttime_ffr)==1 && sign_indicator_ffr==1) ...
        || (sum(ttime(t)==ttime_zlb0)==1 && sign_indicator_zlb0==1) ...
        || (sum(ttime(t)==ttime_zlb1)==1 && sign_indicator_zlb1==1) ...
        || (sum(ttime(t)==ttime_zlb2)==1 && sign_indicator_zlb2==1));

else

    St_old=0;
    DCDhalft_old = eye(nv)*1e8;
end

