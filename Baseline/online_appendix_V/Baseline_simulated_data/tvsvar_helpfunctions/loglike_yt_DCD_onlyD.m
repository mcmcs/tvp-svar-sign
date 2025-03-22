function loglik_prop = loglike_yt_DCD_onlyD(y_t1, B_t1, Q_t1, function_restrictions, ggamma_t1, z_prop)
% likelihood eval y given z_prop
% z_prop = [ddelta, ggamma]

% _onlyD -> z_prop = ddelta


% ---
nv = numel(y_t1);
ddelta_prop = z_prop;
ggamma_prop = ggamma_t1;

% ddelta_prop = z_prop( 1:nv );
% ggamma_prop = z_prop( (nv+1):end );
% MM = size(Vd_random,2);

% ---
%    warning('')
[St_prop, DCDhalft_prop] = function_restrictions(B_t1, ddelta_prop, ggamma_prop, Q_t1);
invSig_half = eye(nv) / DCDhalft_prop';
et = y_t1' * invSig_half;
logyt_prop = -nv/2*log(2*pi) + 1/2*2*sum(log(diag(invSig_half))) -1/2*(et*et');

%    [warnMsg, warnId] = lastwarn;
%     if ~isempty(warnMsg)
% 
%       keyboard
%     end

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

