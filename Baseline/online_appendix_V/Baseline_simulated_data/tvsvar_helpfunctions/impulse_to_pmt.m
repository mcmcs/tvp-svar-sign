function Ih = impulse_to_pmt(Ap_tmp, F1_tmp, F2_tmp, hrz, eqn_ffr, ind_ffr, ind_var)
% -------------------------------------------------------------------------
% impulse response of {response variable} to permanent change in {shock variable}
% Only works for VAR with two lags
% -------------------------------------------------------------------------
% Ap * y(t) = F0 + F1*y(t-1) + F2*y(t-2) + e(t)
% -------------------------------------------------------------------------
% Ap is A'
% L is inv(Ap)
% B is n by k matrix 
% example
% L_tmp  = mat_L0(:,:,tt,gg); %
% Ap_tmp = L_tmp\eye(n); % A'
% B_tmp  = reshape(vec(squeeze(mat_B(tt,:,gg))), [], nvar)';
% B1_tmp = B_tmp(:,2:4); % B1 - coeff on the first lag
% B2_tmp = B_tmp(:,5:7); % B2 - coeff on the second lag
% F1_tmp = (Ap_tmp*B1_tmp); % F1 - coeff on the first lag
% F2_tmp = (Ap_tmp*B2_tmp); % F2 - coeff on the second lag
% -------------------------------------------------------------------------
% compute responses until hrz
% eqn_ffr: location of {response variable} in terms of equation ordering
% ind_ffr: location of {response variable} in terms of variable ordering
% ind_var: location of {shock variable} in terms of variable ordering
% -------------------------------------------------------------------------

% required variables
a_nmz   = Ap_tmp(eqn_ffr, ind_ffr); %normalization
a_shock = -Ap_tmp(eqn_ffr, ind_var)/ a_nmz; % contempareous impact on ffr
f1_ffr  = F1_tmp(eqn_ffr, ind_ffr) / a_nmz; % lagged impact from a ffr
f2_ffr  = F2_tmp(eqn_ffr, ind_ffr) / a_nmz; % lagged impact from a ffr
f1_var  = F1_tmp(eqn_ffr, ind_var) / a_nmz; % lagged impact from a shocking var
f2_var  = F2_tmp(eqn_ffr, ind_var) / a_nmz; % lagged impact from a shocking var

% computing long-run impulses
if ~isinf(hrz)
    Ih = zeros(1,hrz);
    Ih(1) = a_shock;
    Ih(2) = f1_ffr*Ih(1) + f1_var + a_shock;
    Ih(3) = f1_ffr*Ih(2) + f2_ffr*Ih(1) + f1_var + a_shock + f2_var;
    for ih = 4:hrz
        Ih(ih) = f1_ffr*Ih(ih-1) + f2_ffr*Ih(ih-2) + f1_var + a_shock + f2_var;
    end

    % final output (cut Ih in case hrz < 3)
    if hrz <3
        Ih = Ih(1:hrz);
    end
else
    % cumulative effect
    Ih = (f1_var + a_shock + f2_var) / (1-f1_ffr-f2_ffr);
end


