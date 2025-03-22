%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clearvars -global
clear_persistent_variables(fileparts(which('dynare')), false)
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'model_simdata';
M_.dynare_version = '6.1';
oo_.dynare_version = '6.1';
options_.dynare_version = '6.1';
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
M_.endo_names = cell(10,1);
M_.endo_names_tex = cell(10,1);
M_.endo_names_long = cell(10,1);
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
M_.endo_names(8) = {'ytildeobs'};
M_.endo_names_tex(8) = {'ytildeobs'};
M_.endo_names_long(8) = {'ytildeobs'};
M_.endo_names(9) = {'ppiobs'};
M_.endo_names_tex(9) = {'ppiobs'};
M_.endo_names_long(9) = {'ppiobs'};
M_.endo_names(10) = {'Robs'};
M_.endo_names_tex(10) = {'Robs'};
M_.endo_names_long(10) = {'Robs'};
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
M_.endo_nbr = 10;
M_.param_nbr = 14;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.eq_nbr = 10;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
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
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , 'ytilde' ;
  2 , 'name' , 'ppi' ;
  3 , 'name' , 'R' ;
  4 , 'name' , 'rn' ;
  5 , 'name' , 'v' ;
  6 , 'name' , 'a' ;
  7 , 'name' , 'u' ;
  8 , 'name' , 'ytildeobs' ;
  9 , 'name' , 'ppiobs' ;
  10 , 'name' , 'Robs' ;
};
M_.mapping.ytilde.eqidx = [1 2 3 8 ];
M_.mapping.ppi.eqidx = [1 2 3 9 ];
M_.mapping.R.eqidx = [1 3 10 ];
M_.mapping.rn.eqidx = [1 4 ];
M_.mapping.v.eqidx = [3 5 ];
M_.mapping.a.eqidx = [4 6 ];
M_.mapping.u.eqidx = [2 7 ];
M_.mapping.ytildeobs.eqidx = [8 ];
M_.mapping.ppiobs.eqidx = [9 ];
M_.mapping.Robs.eqidx = [10 ];
M_.mapping.vareps_v.eqidx = [5 ];
M_.mapping.vareps_a.eqidx = [6 ];
M_.mapping.vareps_u.eqidx = [7 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 4;
M_.block_structure.block(1).mfs = 4;
M_.block_structure.block(1).equation = [ 5 6 7 4];
M_.block_structure.block(1).variable = [ 5 6 7 4];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 8;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 3 5 6 7 8 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 3;
M_.block_structure.block(2).mfs = 3;
M_.block_structure.block(2).equation = [ 1 2 3];
M_.block_structure.block(2).variable = [ 3 2 1];
M_.block_structure.block(2).is_linear = true;
M_.block_structure.block(2).NNZDerivatives = 11;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [1 4 5 6 8 9 ];
M_.block_structure.block(3).Simulation_Type = 1;
M_.block_structure.block(3).endo_nbr = 3;
M_.block_structure.block(3).mfs = 3;
M_.block_structure.block(3).equation = [ 10 9 8];
M_.block_structure.block(3).variable = [ 10 9 8];
M_.block_structure.block(3).is_linear = true;
M_.block_structure.block(3).NNZDerivatives = 3;
M_.block_structure.block(3).bytecode_jacob_cols_to_sparse = [4 5 6 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([3 1 3 2 3 1 2 3 1 2 1 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 4 4 5 5 6 6 6 8 8 9 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 2 2 2 4 6 9 9 11 12 ]);
M_.block_structure.block(3).g1_sparse_rowval = int32([]);
M_.block_structure.block(3).g1_sparse_colval = int32([]);
M_.block_structure.block(3).g1_sparse_colptr = int32([]);
M_.block_structure.variable_reordered = [ 5 6 7 4 3 2 1 10 9 8];
M_.block_structure.equation_reordered = [ 5 6 7 4 1 2 3 10 9 8];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 3 3;
 5 5;
 6 6;
 7 7;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 1;
 1 3;
 1 4;
 2 1;
 2 2;
 2 7;
 3 1;
 3 2;
 3 3;
 3 5;
 4 4;
 4 6;
 5 5;
 6 6;
 7 7;
 8 1;
 8 8;
 9 2;
 9 9;
 10 3;
 10 10;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 1;
 1 2;
 2 2;
];
M_.block_structure.dyn_tmp_nbr = 0;
M_.state_var = [5 6 7 3 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(14, 1);
M_.endo_trends = struct('deflator', cell(10, 1), 'log_deflator', cell(10, 1), 'growth_factor', cell(10, 1), 'log_growth_factor', cell(10, 1));
M_.NNZDerivatives = [31; -1; -1; ];
M_.dynamic_g1_sparse_rowval = int32([3 5 6 7 1 2 3 8 2 3 9 1 3 10 1 4 3 5 4 6 2 7 8 9 10 1 1 2 5 6 7 ]);
M_.dynamic_g1_sparse_colval = int32([3 5 6 7 11 11 11 11 12 12 12 13 13 13 14 14 15 15 16 16 17 17 18 19 20 21 22 22 31 32 33 ]);
M_.dynamic_g1_sparse_colptr = int32([1 1 1 2 2 3 4 5 5 5 5 9 12 15 17 19 21 23 24 25 26 27 29 29 29 29 29 29 29 29 29 30 31 32 ]);
M_.lhs = {
'ytilde'; 
'ppi'; 
'R'; 
'rn'; 
'v'; 
'a'; 
'u'; 
'ytildeobs'; 
'ppiobs'; 
'Robs'; 
};
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 5];
M_.block_structure_stat.block(1).variable = [ 5];
M_.block_structure_stat.block(2).Simulation_Type = 3;
M_.block_structure_stat.block(2).endo_nbr = 1;
M_.block_structure_stat.block(2).mfs = 1;
M_.block_structure_stat.block(2).equation = [ 6];
M_.block_structure_stat.block(2).variable = [ 6];
M_.block_structure_stat.block(3).Simulation_Type = 3;
M_.block_structure_stat.block(3).endo_nbr = 1;
M_.block_structure_stat.block(3).mfs = 1;
M_.block_structure_stat.block(3).equation = [ 7];
M_.block_structure_stat.block(3).variable = [ 7];
M_.block_structure_stat.block(4).Simulation_Type = 1;
M_.block_structure_stat.block(4).endo_nbr = 1;
M_.block_structure_stat.block(4).mfs = 1;
M_.block_structure_stat.block(4).equation = [ 4];
M_.block_structure_stat.block(4).variable = [ 4];
M_.block_structure_stat.block(5).Simulation_Type = 6;
M_.block_structure_stat.block(5).endo_nbr = 3;
M_.block_structure_stat.block(5).mfs = 3;
M_.block_structure_stat.block(5).equation = [ 1 2 3];
M_.block_structure_stat.block(5).variable = [ 3 1 2];
M_.block_structure_stat.block(6).Simulation_Type = 1;
M_.block_structure_stat.block(6).endo_nbr = 3;
M_.block_structure_stat.block(6).mfs = 3;
M_.block_structure_stat.block(6).equation = [ 10 9 8];
M_.block_structure_stat.block(6).variable = [ 10 9 8];
M_.block_structure_stat.variable_reordered = [ 5 6 7 4 3 1 2 10 9 8];
M_.block_structure_stat.equation_reordered = [ 5 6 7 4 1 2 3 10 9 8];
M_.block_structure_stat.incidence.sparse_IM = [
 1 2;
 1 3;
 1 4;
 2 1;
 2 2;
 2 7;
 3 1;
 3 2;
 3 3;
 3 5;
 4 4;
 4 6;
 5 5;
 6 6;
 7 7;
 8 1;
 8 8;
 9 2;
 9 9;
 10 3;
 10 10;
];
M_.block_structure_stat.tmp_nbr = 0;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(4).g1_sparse_rowval = int32([]);
M_.block_structure_stat.block(4).g1_sparse_colval = int32([]);
M_.block_structure_stat.block(4).g1_sparse_colptr = int32([]);
M_.block_structure_stat.block(5).g1_sparse_rowval = int32([1 3 2 3 1 2 3 ]);
M_.block_structure_stat.block(5).g1_sparse_colval = int32([1 1 2 2 3 3 3 ]);
M_.block_structure_stat.block(5).g1_sparse_colptr = int32([1 3 5 8 ]);
M_.block_structure_stat.block(6).g1_sparse_rowval = int32([]);
M_.block_structure_stat.block(6).g1_sparse_colval = int32([]);
M_.block_structure_stat.block(6).g1_sparse_colptr = int32([]);
M_.static_g1_sparse_rowval = int32([2 3 8 1 2 3 9 1 3 10 1 4 3 5 4 6 2 7 8 9 10 ]);
M_.static_g1_sparse_colval = int32([1 1 1 2 2 2 2 3 3 3 4 4 5 5 6 6 7 7 8 9 10 ]);
M_.static_g1_sparse_colptr = int32([1 4 8 11 13 15 17 19 20 21 22 ]);
M_.params(1) = 0.99;
bbeta = M_.params(1);
M_.params(2) = 0.10;
kkappa = M_.params(2);
M_.params(3) = 1.00;
inv_ssigma = M_.params(3);
M_.params(4) = 0.85;
rrho_i = M_.params(4);
M_.params(5) = 1.50;
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
oo_.steady_state(8) = 0;
oo_.steady_state(9) = 0;
oo_.steady_state(10) = 400*oo_.steady_state(3);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
options_resid_ = struct();
display_static_residuals(M_, options_, oo_, options_resid_);
oo_.dr.eigval = check(M_,options_,oo_);
steady;
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
M_.Sigma_e(3, 3) = (1)^2;
options_.drop = 0;
options_.irf = 0;
options_.nograph = true;
options_.noprint = true;
options_.order = 1;
options_.periods = 215;
var_list_ = {'ytildeobs';'ppiobs';'Robs'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model_simdata_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
