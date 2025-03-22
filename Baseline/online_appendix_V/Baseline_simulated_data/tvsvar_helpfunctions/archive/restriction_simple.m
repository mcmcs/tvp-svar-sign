function [St_old, DCDhalft_old] = restriction_simple(B_t, ddelta_t, ggamma_t, Q_t,extra_arg)

% extra_arg: for the backward compatibility

% simple sign restriction function
% this is for DGP to do Geweke test
Dt_old    = diag(exp(ddelta_t/2));
Ct_old    = ggammatoC(ggamma_t);
DCDhalft_old = chol(Dt_old*Ct_old*Dt_old, 'lower');
Lt_old = DCDhalft_old * Q_t;
St_old = (Lt_old(1,1) > 0) & (Lt_old(2,1) < 0);