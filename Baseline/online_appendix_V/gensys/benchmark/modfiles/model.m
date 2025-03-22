%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'model';
M_.dynare_version = '4.5.7';
oo_.dynare_version = '4.5.7';
options_.dynare_version = '4.5.7';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('model.log');
M_.exo_names = 'vareps_v';
M_.exo_names_tex = 'vareps\_v';
M_.exo_names_long = 'vareps_v';
M_.exo_names = char(M_.exo_names, 'vareps_a');
M_.exo_names_tex = char(M_.exo_names_tex, 'vareps\_a');
M_.exo_names_long = char(M_.exo_names_long, 'vareps_a');
M_.exo_names = char(M_.exo_names, 'vareps_u');
M_.exo_names_tex = char(M_.exo_names_tex, 'vareps\_u');
M_.exo_names_long = char(M_.exo_names_long, 'vareps_u');
M_.endo_names = 'ytilde';
M_.endo_names_tex = 'ytilde';
M_.endo_names_long = 'ytilde';
M_.endo_names = char(M_.endo_names, 'ppi');
M_.endo_names_tex = char(M_.endo_names_tex, 'ppi');
M_.endo_names_long = char(M_.endo_names_long, 'ppi');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'rn');
M_.endo_names_tex = char(M_.endo_names_tex, 'rn');
M_.endo_names_long = char(M_.endo_names_long, 'rn');
M_.endo_names = char(M_.endo_names, 'v');
M_.endo_names_tex = char(M_.endo_names_tex, 'v');
M_.endo_names_long = char(M_.endo_names_long, 'v');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'u');
M_.endo_names_tex = char(M_.endo_names_tex, 'u');
M_.endo_names_long = char(M_.endo_names_long, 'u');
M_.endo_names = char(M_.endo_names, 'ytildeobs');
M_.endo_names_tex = char(M_.endo_names_tex, 'ytildeobs');
M_.endo_names_long = char(M_.endo_names_long, 'ytildeobs');
M_.endo_names = char(M_.endo_names, 'ppiobs');
M_.endo_names_tex = char(M_.endo_names_tex, 'ppiobs');
M_.endo_names_long = char(M_.endo_names_long, 'ppiobs');
M_.endo_names = char(M_.endo_names, 'Robs');
M_.endo_names_tex = char(M_.endo_names_tex, 'Robs');
M_.endo_names_long = char(M_.endo_names_long, 'Robs');
M_.endo_partitions = struct();
M_.param_names = 'bbeta';
M_.param_names_tex = 'bbeta';
M_.param_names_long = 'bbeta';
M_.param_names = char(M_.param_names, 'kkappa');
M_.param_names_tex = char(M_.param_names_tex, 'kkappa');
M_.param_names_long = char(M_.param_names_long, 'kkappa');
M_.param_names = char(M_.param_names, 'inv_ssigma');
M_.param_names_tex = char(M_.param_names_tex, 'inv\_ssigma');
M_.param_names_long = char(M_.param_names_long, 'inv_ssigma');
M_.param_names = char(M_.param_names, 'rrho_i');
M_.param_names_tex = char(M_.param_names_tex, 'rrho\_i');
M_.param_names_long = char(M_.param_names_long, 'rrho_i');
M_.param_names = char(M_.param_names, 'pphi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'pphi\_pi');
M_.param_names_long = char(M_.param_names_long, 'pphi_pi');
M_.param_names = char(M_.param_names, 'pphi_ytilde');
M_.param_names_tex = char(M_.param_names_tex, 'pphi\_ytilde');
M_.param_names_long = char(M_.param_names_long, 'pphi_ytilde');
M_.param_names = char(M_.param_names, 'ppsi_yan');
M_.param_names_tex = char(M_.param_names_tex, 'ppsi\_yan');
M_.param_names_long = char(M_.param_names_long, 'ppsi_yan');
M_.param_names = char(M_.param_names, 'rrho_v');
M_.param_names_tex = char(M_.param_names_tex, 'rrho\_v');
M_.param_names_long = char(M_.param_names_long, 'rrho_v');
M_.param_names = char(M_.param_names, 'rrho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rrho\_a');
M_.param_names_long = char(M_.param_names_long, 'rrho_a');
M_.param_names = char(M_.param_names, 'rrho_u');
M_.param_names_tex = char(M_.param_names_tex, 'rrho\_u');
M_.param_names_long = char(M_.param_names_long, 'rrho_u');
M_.param_names = char(M_.param_names, 'ssigma_v');
M_.param_names_tex = char(M_.param_names_tex, 'ssigma\_v');
M_.param_names_long = char(M_.param_names_long, 'ssigma_v');
M_.param_names = char(M_.param_names, 'ssigma_a');
M_.param_names_tex = char(M_.param_names_tex, 'ssigma\_a');
M_.param_names_long = char(M_.param_names_long, 'ssigma_a');
M_.param_names = char(M_.param_names, 'ssigma_u');
M_.param_names_tex = char(M_.param_names_tex, 'ssigma\_u');
M_.param_names_long = char(M_.param_names_long, 'ssigma_u');
M_.param_names = char(M_.param_names, 'rrho');
M_.param_names_tex = char(M_.param_names_tex, 'rrho');
M_.param_names_long = char(M_.param_names_long, 'rrho');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 10;
M_.param_nbr = 14;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'ytildeobs'};
options_.varobs(2)  = {'Robs'};
options_.varobs(3)  = {'ppiobs'};
options_.varobs_id = [ 8 10 9  ];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('model_static');
erase_compiled_function('model_dynamic');
M_.orig_eq_nbr = 10;
M_.eq_nbr = 10;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 5 15;
 0 6 16;
 1 7 0;
 0 8 0;
 2 9 0;
 3 10 0;
 4 11 0;
 0 12 0;
 0 13 0;
 0 14 0;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 4;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(14, 1);
M_.NNZDerivatives = [31; 0; -1];
M_.params( 1 ) = 0.99;
bbeta = M_.params( 1 );
M_.params( 2 ) = 0.25;
kkappa = M_.params( 2 );
M_.params( 3 ) = 2.00;
inv_ssigma = M_.params( 3 );
M_.params( 4 ) = 0.5;
rrho_i = M_.params( 4 );
M_.params( 5 ) = 2.0;
pphi_pi = M_.params( 5 );
M_.params( 6 ) = 0.25;
pphi_ytilde = M_.params( 6 );
M_.params( 7 ) = 1;
ppsi_yan = M_.params( 7 );
M_.params( 8 ) = 0.5;
rrho_v = M_.params( 8 );
M_.params( 9 ) = 0.5;
rrho_a = M_.params( 9 );
M_.params( 10 ) = 0.5;
rrho_u = M_.params( 10 );
switch datasets{i_var}
case 'baseline'
set_param_value('ssigma_v',0.25/100);
set_param_value('ssigma_a',5/100);
set_param_value('ssigma_u',0.25/100);
end
M_.params( 14 ) = (-log(M_.params(1)));
rrho = M_.params( 14 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 0;
oo_.steady_state( 2 ) = 0;
oo_.steady_state( 3 ) = M_.params(14);
oo_.steady_state( 4 ) = M_.params(14);
oo_.steady_state( 5 ) = 0;
oo_.steady_state( 6 ) = 0;
oo_.steady_state( 7 ) = 0;
oo_.steady_state( 8 ) = 0;
oo_.steady_state( 9 ) = 0;
oo_.steady_state( 10 ) = 400*oo_.steady_state(3);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid;
oo_.dr.eigval = check(M_,options_,oo_);
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
M_.Sigma_e(3, 3) = (1)^2;
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 2, 0.25, 1e-5, 10, 2, 0.25, 0.20, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 3, 2.0, 1e-5, 10, 2, 2, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, 0.5, 1e-5, 1, 1, 0.5000, 0.2000, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 5, 2.0, 1e-5, 10, 2, 2.0, 0.50, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, 0.25, 1e-5, 10, 2, 0.25, 0.20, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, 0.5, 1e-5, 1, 1, 0.50, 0.20, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, 0.5, 1e-5, 1, 1, 0.50, 0.20, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 10, 0.5, 1e-5, 1, 1, 0.50, 0.20, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 11, 0.01, 1e-5, 10, 6, 0.01, 0.50, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 12, 0.01, 1e-5, 10, 6, 0.01, 0.50, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 13, 0.01, 1e-5, 10, 6, 0.01, 0.50, NaN, NaN, NaN ];
resid;
oo_.dr.eigval = check(M_,options_,oo_);
options_.filtered_vars = 1;
options_.mh_drop = 0.25;
options_.mh_jscale = 0.6;
options_.mh_nblck = 1;
options_.mh_replic = 50000;
options_.mode_compute = 6;
options_.plot_priors = 0;
options_.prior_trunc = 0;
options_.MCMC_jumping_covariance = 'hessian';
options_.datafile = '~/Dropbox/nk3/estimatedNK3_public/simdata.mat';
options_.endo_vars_for_moment_computations_in_estimation = 'all_endogenous_variables';
options_.filter_step_ahead = 1;
options_.first_obs = 1;
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
save('model_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('model_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('model_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('model_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('model_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('model_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('model_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
