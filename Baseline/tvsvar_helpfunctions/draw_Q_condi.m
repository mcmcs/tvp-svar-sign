function Q_tmp = draw_Q_condi(B_t, ddelta_t, ggamma_t, function_restrictions, info_nvar, t)
% draw Q from conditional posterior
St_old = false;
while (~St_old)
    Q_tmp = draw_from_On(info_nvar);
    St_old = function_restrictions(B_t, ddelta_t, ggamma_t, Q_tmp, t);
end
