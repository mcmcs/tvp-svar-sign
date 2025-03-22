function loglik_prop = loglike_yt_B(y_t1, Z_t1, ddelta_t1, ggamma_t1, Q_t1, function_restrictions, z_prop)
% likelihood eval y given z_prop
% z_prop = bbet


% ---
nv = numel(y_t1);
bbet_prop = z_prop;

% ---
[St_prop, DCDhalft_prop] = function_restrictions(bbet_prop, ddelta_t1, ggamma_t1, Q_t1);
invSig_half = eye(nv) / DCDhalft_prop';
et = (y_t1 - Z_t1*bbet_prop)' * invSig_half;
logyt_prop = -nv/2*log(2*pi) + 1/2*2*sum(log(diag(invSig_half))) -1/2*(et*et');

% For joint prior we do not need to compute normalizing constant
logRt_prop = 0;

% --- normalizing constant part
% NOTE: we don't have normalizing constant in the last period
% if ~is_last
%     ddelta_prop2 = ddelta_prop + Vd_half * Vd_random;
%     ggamma_prop2 = ggamma_prop + Vg_half * Vg_random;
%     Q_prop2 = Q_random;
%     Rt_j = zeros(MM,1);
%     for j2 = 1:MM
%         St2_old = function_restrictions(ddelta_prop2(:,j2), ggamma_prop2(:,j2), Q_prop2(:,:,j2));
%         Rt_j(j2) = St2_old;
%     end
%     logRt_prop = log( mean(Rt_j) );
% else
%     logRt_prop =0;
% end

loglik_prop = -Inf;
if St_prop
    loglik_prop = logyt_prop - logRt_prop;
end

