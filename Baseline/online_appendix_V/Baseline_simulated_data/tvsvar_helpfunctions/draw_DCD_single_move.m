function [ddelta_old, ggamma_old, ntry_old] = draw_DCD_single_move(yhat_old, B_old, ddelta_old, ggamma_old, Q_old, Vd_old, Vg_old, ddelta_hat, ddelta_hat_Vinv, ggamma_hat, ggamma_hat_Vinv, function_restrictions)

% function to implement "single-move" sampling together with elliptical
% slice sampling

% we include B_old as an argument because function_restrictions may depend
% on B_old

[T,nv] = size(yhat_old);
ng = nv * (nv-1)/2;
ntry_old = zeros(T,1);

% --- Other parameters
% Vd      = diag(Vd_old);
% Vd_half = chol(Vd,'lower');
% Vg      = diag(Vg_old);
% Vg_half = chol(diag(Vg_old),'lower');

Vd      = diag(Vd_old);
% Vd_half = diag(sqrt(Vd_old));
Vg      = diag(Vg_old);
% Vg_half = diag(sqrt(Vg_old));

% --- Loop over time
for t=1:T
    
    % build relevant data at time t
    ddelta_t1 = ddelta_old(t,:)';
    ggamma_t1 = ggamma_old(t,:)';
    
    
    if t < T
        ddelta_t2 = ddelta_old(t+1,:)';
        ggamma_t2 = ggamma_old(t+1,:)';
    else
        ddelta_t2 = [];
        ggamma_t2 = [];
    end
    
    if t == 1 % ddelta_t0 and ggamma_t0 are only used to calculate prior
        ddelta_t0 = ddelta_hat;
        ggamma_t0 = ggamma_hat;
    else
        ddelta_t0 = ddelta_old(t-1,:)';
        ggamma_t0 = ggamma_old(t-1,:)';
    end
    
    % --- Prior
    
    % prior mean = (1/2) * (gamma(t+1) + gamma(t))
    % prior variance = (1/2) * Vg
    
    if t == T
        ggamma_pr_m = ggamma_t0;
        ggamma_pr_V = Vg;
        ddelta_pr_m = ddelta_t0;
        ddelta_pr_V = Vd;
        is_last = true;
    elseif t == 1
        Vg_inv = diag(Vg_old.^-1); %Vg_old = diag(Vg)
        ggamma_pr_V = eye(ng) / ( ggamma_hat_Vinv +  Vg_inv ) ;
        ggamma_pr_m = ggamma_pr_V * ( Vg_inv*ggamma_t2 + ggamma_hat_Vinv*ggamma_t0);
        
        Vd_inv = diag(Vd_old.^-1); %Vd_old = diag(Vd)
        ddelta_pr_V = eye(nv) / ( ddelta_hat_Vinv +  Vd_inv ) ;
        ddelta_pr_m = ddelta_pr_V * ( Vd_inv*ddelta_t2 + ddelta_hat_Vinv*ddelta_t0);
        is_last = false;
    else
        ggamma_pr_m = (1/2) * ( ggamma_t2 + ggamma_t0);
        ggamma_pr_V = (1/2) * Vg;
        ddelta_pr_m = (1/2) * ( ddelta_t2 + ddelta_t0);
        ddelta_pr_V = (1/2) * Vd;
        is_last = false;
    end
    % --- needs to lik eval
    y_t1 = yhat_old(t,:)';
    Q_t1 = Q_old(:,:,t);
    B_t1 = B_old(t,:);
    
    % --- for slice sampling
    function_restrictions_t = @(ag0, ag1, ag2, ag3) function_restrictions(ag0, ag1, ag2, ag3, t); %time-varying restrictions
    
    slice.scale_z    = 1;
    slice.mean       = [ddelta_pr_m; ggamma_pr_m];
    slice.chol_cov_z = chol(blkdiag(ddelta_pr_V, ggamma_pr_V),'lower');
    slice.fcn_lik    = @(z_prop) loglike_yt_DCD(y_t1, B_t1, Q_t1, function_restrictions_t, z_prop);
    slice.nobs       = nv + ng;
    
    % --- actual slice sampling
    z_old = [ddelta_t1; ggamma_t1];
    lik_old = slice.fcn_lik(z_old);
    
    [z_old, ~, n_try] = slice_sampling_v02(slice, z_old, lik_old);
    
    % --- unroll
    ddelta_old(t,:) = z_old(1:nv)';
    ggamma_old(t,:) = z_old((nv+1):end)';
    ntry_old(t,:) = n_try;
    
end





