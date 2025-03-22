function [model_ss,model_para_ss] = get_model_ss(ttheta)


calibrated_parameters = get_calibrated_parameters();

%% calibrated parameters
bbeta = calibrated_parameters(1,1);
%inv_ssigma =calibrated_parameters(2,1);
%rrho_i = calibrated_parameters(3,1);
%pphi_pi = calibrated_parameters(4,1);
%pphi_ytilde = calibrated_parameters(5,1);
%ppsi_yan = calibrated_parameters(6,1);
%rrho_v = calibrated_parameters(7,1);
%rrho_a= calibrated_parameters(8,1);
%rrho_u = calibrated_parameters(9,1);
%ssigma_v = calibrated_parameters(10,1);
%ssigma_a = calibrated_parameters(11,1);
%ssigma_u =calibrated_parameters(12,1);

%% estimated/para of interest parameters

% Structural parameters
% kkappa          = ttheta(1);




%% steady state and matrices for GENSYS




rrho         = -log(bbeta);






ytilde    = 0;   
ppi       = 0;   
R         = rrho;  
rn        = rrho;  
v         = 0;    
a         = 0;    
u         = 0;      

model_ss(1) = ytilde;
model_ss(2) = ppi;
model_ss(3) = R;         
model_ss(4) = rn;
model_ss(5) = v;
model_ss(6) = a;
model_ss(7) = u;



model_para_ss(1,1) = rrho;



end

