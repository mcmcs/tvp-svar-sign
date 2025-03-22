function [B_old, ntry_old] = draw_B_single_move3(y_old, x_old, B_old, ddelta_old, ggamma_old, Q_old, VB_old, bbet_hat, bbet_hat_Vinv, function_restrictions)

% move3: elliptical sampling

% function to implement "single-move" sampling together with elliptical
% slice sampling

[T,nv] = size(y_old);
ng = nv * (nv-1)/2;
ntry_old = zeros(T,1);
k = size(VB_old,1);
eye_k = eye(k);

% --- Other parameters


% --- Loop over time
for t=1:T
    
    % build relevant data at time t
    bbet_t1 = B_old(t,:)';
    
    if t < T
        bbet_t2 = B_old(t+1,:)';
    else
        bbet_t2 = [];
    end
    
    if t == 1 % ddelta_t0 and ggamma_t0 are only used to calculate prior
        bbet_t0 = bbet_hat';
    else
        bbet_t0 = B_old(t-1,:)';
    end
    
    % --- Prior
    
    % prior mean = (1/2) * (gamma(t+1) + gamma(t))
    % prior variance = (1/2) * Vg
    
    if t == T
        bbet_pr_m = bbet_t0;
        bbet_pr_V = VB_old;
    elseif t == 1
        VB_inv = eye(k) / VB_old;
        bbet_pr_V = eye(k) / ( bbet_hat_Vinv + VB_inv);
        
%         VB_inv = inv(VB_old);
%         bbet_pr_V = inv( bbet_hat_Vinv + VB_inv);
        
        bbet_pr_m = bbet_pr_V * (VB_inv*bbet_t2 + bbet_hat_Vinv*bbet_t0);
        
%         bbet_hat_V = inv(bbet_hat_Vinv);
%         
%         kf_sigma = bbet_hat_V + VB_old; %Z_t1*bbet_pr_V*Z_t1' + DCD_t1;
%         kf_k = bbet_hat_V / kf_sigma; %(bbet_pr_V*Z_t1')/kf_sigma;
%         kf_ky = kf_k * bbet_t2; %kf_k *y_t1;
%         kf_IkZ = eye_k - kf_k;
% 
%         bbet_pr_m = kf_IkZ * bbet_t0 + kf_ky;
%         bbet_pr_V = bbet_hat_V - kf_k * bbet_hat_V;

    else
        bbet_pr_m = (1/2) * ( bbet_t2 + bbet_t0);
        bbet_pr_V = (1/2) * VB_old;
    end
    % --- needs to lik eval
    y_t1 = y_old(t,:)';
    %Z_t1 = Z_old([0:nv-1]*T+t,:);
    Z_t1 = kron( eye(nv), x_old(t,:) );
    Q_t1 = Q_old(:,:,t);
    ddelta_t1 = ddelta_old(t,:);
    ggamma_t1 = ggamma_old(t,:);
%     Dt_old    = diag(exp(ddelta_old(t,:)/2));
%     Ct_old    = ggammatoC(ggamma_old(t,:));
%     DCD_t1 = Dt_old*Ct_old*Dt_old;

    % --- for slice sampling
    function_restrictions_t = @(ag0, ag1, ag2, ag3) function_restrictions(ag0, ag1, ag2, ag3, t); %time-varying restrictions
    
    slice.scale_z = 1;
    slice.mean = bbet_pr_m;
    slice.chol_cov_z = chol(bbet_pr_V,'lower');
    slice.fcn_lik = @(z_prop) loglike_yt_B(y_t1, Z_t1, ddelta_t1, ggamma_t1, Q_t1, function_restrictions_t, z_prop);
    slice.nobs = k;
    
    % --- actual slice sampling
    z_old = bbet_t1;
    lik_old = slice.fcn_lik(z_old);
    [bbet_t1_old, ~, n_try] = slice_sampling_v02(slice, z_old, lik_old);
    
    
	% --- unroll
    B_old(t,:) = bbet_t1_old;
    ntry_old(t,:) = n_try;
    
    
    % --- Post
%     kf_sigma = Z_t1*bbet_pr_V*Z_t1' + DCD_t1;
%     kf_k     = (bbet_pr_V*Z_t1')/kf_sigma;
%     kf_ky    = kf_k *y_t1;
%     kf_IkZ   = eye_k - kf_k*Z_t1;
% 
%     bbet_post_m = kf_IkZ * bbet_pr_m + kf_ky; % posterior mean
%     bbet_post_V = bbet_pr_V - kf_k * Z_t1 * bbet_pr_V;
%     bbet_post_Vhalf = chol(bbet_post_V, 'lower');
    
    
%     kf_sigma = Z_t1*bbet_pr_V*Z_t1' + DCD_t1;
%     kf_k = (bbet_pr_V*Z_t1')/kf_sigma;
%     kf_ky = kf_k *y_t1;
%     kf_IkZ = eye_k - kf_k*Z_t1;
%     
%     bbet_post_m = kf_IkZ * bbet_pr_m + kf_ky; %only this depends on particles
%     bbet_post_V = bbet_pr_V - kf_k * Z_t1*bbet_pr_V;
%     bbet_post_Vhalf = chol(bbet_post_V,'lower');
            
%     ntry = 0;
%     St_old = false;
%     while (~St_old)
%         bbet_t1 = ( bbet_post_m + bbet_post_Vhalf * randn(k,1) )';
%         St_old = function_restrictions(bbet_t1, ddelta_t1, ggamma_t1, Q_t1, t);
%         ntry = ntry + 1;
%         
%         
%         if ntry > 1000
%             keyboard;
%         end
%     end
    
%     % --- unroll
%     B_old(t,:) = bbet_t1;
%     ntry_old(t,:) = ntry;

    
end





