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
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'model';
M_.dynare_version = '5.0';
oo_.dynare_version = '5.0';
options_.dynare_version = '5.0';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'vareps_v'};
M_.exo_names_tex(1) = {'vareps\_v'};
M_.exo_names_long(1) = {'vareps_v'};
M_.exo_names(2) = {'vareps_a'};
M_.exo_names_tex(2) = {'vareps\_a'};
M_.exo_names_long(2) = {'vareps_a'};
M_.exo_names(3) = {'vareps_u'};
M_.exo_names_tex(3) = {'vareps\_u'};
M_.exo_names_long(3) = {'vareps_u'};
M_.endo_names = cell(7,1);
M_.endo_names_tex = cell(7,1);
M_.endo_names_long = cell(7,1);
M_.endo_names(1) = {'ytilde'};
M_.endo_names_tex(1) = {'ytilde'};
M_.endo_names_long(1) = {'ytilde'};
M_.endo_names(2) = {'ppi'};
M_.endo_names_tex(2) = {'ppi'};
M_.endo_names_long(2) = {'ppi'};
M_.endo_names(3) = {'R'};
M_.endo_names_tex(3) = {'R'};
M_.endo_names_long(3) = {'R'};
M_.endo_names(4) = {'rn'};
M_.endo_names_tex(4) = {'rn'};
M_.endo_names_long(4) = {'rn'};
M_.endo_names(5) = {'v'};
M_.endo_names_tex(5) = {'v'};
M_.endo_names_long(5) = {'v'};
M_.endo_names(6) = {'a'};
M_.endo_names_tex(6) = {'a'};
M_.endo_names_long(6) = {'a'};
M_.endo_names(7) = {'u'};
M_.endo_names_tex(7) = {'u'};
M_.endo_names_long(7) = {'u'};
M_.endo_partitions = struct();
M_.param_names = cell(14,1);
M_.param_names_tex = cell(14,1);
M_.param_names_long = cell(14,1);
M_.param_names(1) = {'bbeta'};
M_.param_names_tex(1) = {'bbeta'};
M_.param_names_long(1) = {'bbeta'};
M_.param_names(2) = {'kkappa'};
M_.param_names_tex(2) = {'kkappa'};
M_.param_names_long(2) = {'kkappa'};
M_.param_names(3) = {'inv_ssigma'};
M_.param_names_tex(3) = {'inv\_ssigma'};
M_.param_names_long(3) = {'inv_ssigma'};
M_.param_names(4) = {'rrho_i'};
M_.param_names_tex(4) = {'rrho\_i'};
M_.param_names_long(4) = {'rrho_i'};
M_.param_names(5) = {'pphi_pi'};
M_.param_names_tex(5) = {'pphi\_pi'};
M_.param_names_long(5) = {'pphi_pi'};
M_.param_names(6) = {'pphi_ytilde'};
M_.param_names_tex(6) = {'pphi\_ytilde'};
M_.param_names_long(6) = {'pphi_ytilde'};
M_.param_names(7) = {'ppsi_yan'};
M_.param_names_tex(7) = {'ppsi\_yan'};
M_.param_names_long(7) = {'ppsi_yan'};
M_.param_names(8) = {'rrho_v'};
M_.param_names_tex(8) = {'rrho\_v'};
M_.param_names_long(8) = {'rrho_v'};
M_.param_names(9) = {'rrho_a'};
M_.param_names_tex(9) = {'rrho\_a'};
M_.param_names_long(9) = {'rrho_a'};
M_.param_names(10) = {'rrho_u'};
M_.param_names_tex(10) = {'rrho\_u'};
M_.param_names_long(10) = {'rrho_u'};
M_.param_names(11) = {'ssigma_v'};
M_.param_names_tex(11) = {'ssigma\_v'};
M_.param_names_long(11) = {'ssigma_v'};
M_.param_names(12) = {'ssigma_a'};
M_.param_names_tex(12) = {'ssigma\_a'};
M_.param_names_long(12) = {'ssigma_a'};
M_.param_names(13) = {'ssigma_u'};
M_.param_names_tex(13) = {'ssigma\_u'};
M_.param_names_long(13) = {'ssigma_u'};
M_.param_names(14) = {'rrho'};
M_.param_names_tex(14) = {'rrho'};
M_.param_names_long(14) = {'rrho'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 7;
M_.param_nbr = 14;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 7;
M_.eq_nbr = 7;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 5 12;
 0 6 13;
 1 7 0;
 0 8 0;
 2 9 0;
 3 10 0;
 4 11 0;]';
M_.nstatic = 1;
M_.nfwrd   = 2;
M_.npred   = 4;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 4;
M_.ndynamic   = 6;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'ytilde' ;
  2 , 'name' , 'ppi' ;
  3 , 'name' , 'R' ;
  4 , 'name' , 'rn' ;
  5 , 'name' , 'v' ;
  6 , 'name' , 'a' ;
  7 , 'name' , 'u' ;
};
M_.mapping.ytilde.eqidx = [1 2 3 ];
M_.mapping.ppi.eqidx = [1 2 3 ];
M_.mapping.R.eqidx = [1 3 ];
M_.mapping.rn.eqidx = [1 4 ];
M_.mapping.v.eqidx = [3 5 ];
M_.mapping.a.eqidx = [4 6 ];
M_.mapping.u.eqidx = [2 7 ];
M_.mapping.vareps_v.eqidx = [5 ];
M_.mapping.vareps_a.eqidx = [6 ];
M_.mapping.vareps_u.eqidx = [7 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 5 6 7 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(14, 1);
M_.endo_trends = struct('deflator', cell(7, 1), 'log_deflator', cell(7, 1), 'growth_factor', cell(7, 1), 'log_growth_factor', cell(7, 1));
M_.NNZDerivatives = [25; -1; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.99;
bbeta = M_.params(1);
M_.params(2) = 0.5;
kkappa = M_.params(2);
M_.params(3) = 1.00;
inv_ssigma = M_.params(3);
M_.params(4) = 0.85;
rrho_i = M_.params(4);
M_.params(5) = 1.5;
pphi_pi = M_.params(5);
M_.params(6) = 0.125;
pphi_ytilde = M_.params(6);
M_.params(7) = 1;
ppsi_yan = M_.params(7);
M_.params(8) = 0.0;
rrho_v = M_.params(8);
M_.params(9) = 0.9;
rrho_a = M_.params(9);
M_.params(10) = 0.9;
rrho_u = M_.params(10);
set_param_value('ssigma_v',0.1/100);
set_param_value('ssigma_a',0.05/100);
set_param_value('ssigma_u',0.05/100);
M_.params(14) = (-log(M_.params(1)));
rrho = M_.params(14);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = 0;
oo_.steady_state(2) = 0;
oo_.steady_state(3) = M_.params(14);
oo_.steady_state(4) = M_.params(14);
oo_.steady_state(5) = 0;
oo_.steady_state(6) = 0;
oo_.steady_state(7) = 0;
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
resid;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 20;
options_.order = 1;
var_list_ = {'ytilde';'ppi';'R'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
switch kkappa
case 0.1
oo_flat = oo_;
save('flat_PC.mat','oo_flat')
case 0.5
oo_steep = oo_;
save('steep_PC.mat','oo_steep')
end


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
