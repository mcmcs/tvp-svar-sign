function [out,errors,model_ss,dynare_ss] = get_dynare_check(tol,ttheta0,oo_)
% check that the model solution match dynare


[GAM0, GAM1, PPSI,PPI,C] = get_model(ttheta0);                      % Linearized model evaluated at ttheta0
[T1,~,T0,~,~,~,~,~,~] = gensys(GAM0,GAM1,C,PPSI,PPI); % Policy functions evaluated at ttheta0


% Robustness check: compare policy functions with dynare solutions:
% Compare steady state
dynare_ss = oo_.steady_state;
[model_ss,~] = get_model_ss(ttheta0);
model_ss  = model_ss';
err_ss    = max(abs(model_ss - dynare_ss));

% Compare policy functions
n                      = size(T1,1);
dynare_order           = oo_.dr.inv_order_var';
dynare_displayed_T1    = oo_.dr.ghx(dynare_order,:)';
dynare_displayed_T0    = oo_.dr.ghu(dynare_order,:)';
dynare_state_var       = oo_.dr.state_var';

% Determine the variables in Gensys state vector for which Dynare does not
% compute policy functions, i.e. the terms E_{t}(x_{t+1}) for any x.
gensys_forward_looking    = (8:9); % input by user: model specific
dynare_exo_var            = [5 6 7]'; % input by user: model specific

dynare_pf_to_gensys_T1  = zeros(n,n);
dynare_pf_to_gensys_T0  = zeros(n,size(T0,2));

for j=1:size(dynare_state_var,1)
    if dynare_state_var(j,1)<=(gensys_forward_looking(1)-1)
    dynare_pf_to_gensys_T1(:,dynare_state_var(j,1))    = [dynare_displayed_T1(j,1:(gensys_forward_looking(1)-1)),zeros(1,size(gensys_forward_looking ,2)),dynare_displayed_T1(j,gensys_forward_looking(1):end)]';
    else
    dynare_pf_to_gensys_T1(:,11+dynare_state_var(j,1)) = [dynare_displayed_T1(j,1:(gensys_forward_looking(1)-1)),zeros(1,size(gensys_forward_looking ,2)),dynare_displayed_T1(j,gensys_forward_looking(1):end)]';    
    end
end
dynare_pf_to_gensys_T1(gensys_forward_looking,:)    = T1(gensys_forward_looking,:);
diff_T1                                             = T1- dynare_pf_to_gensys_T1;
err_T1                                              = max(max(abs(diff_T1)));
for j=1:size(dynare_exo_var,1)
    dynare_pf_to_gensys_T0(:,j) = [dynare_displayed_T0(j,1:(gensys_forward_looking(1)-1)),zeros(1,size(gensys_forward_looking ,2)),dynare_displayed_T0(j,gensys_forward_looking(1):end)]';
end
dynare_pf_to_gensys_T0(gensys_forward_looking,:)    = T0(gensys_forward_looking,:);
diff_T0 = T0 - dynare_pf_to_gensys_T0;
err_T0  = max(max(abs(diff_T0)));

errors = [err_ss,err_T1,err_T0];

if max(max(max(errors)))<tol
    out =1;
else
    out = 0; 
end




end