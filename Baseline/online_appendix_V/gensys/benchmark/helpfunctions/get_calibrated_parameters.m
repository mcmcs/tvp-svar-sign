function cal_para = get_calibrated_parameters()

%% calibrated parameters

bbeta        = 0.99; % discount factor
kkappa       = 0.10; % slope of the Phillips curve
inv_ssigma   = 1.00; % elasticity of intertemporal substitution
rrho_i       = 0.85; % monetary policy inertia
pphi_pi      = 1.50; % response to inflation
pphi_ytilde  = 0.5/4; % response to the output gap

 
   
    

ppsi_yan     = 1; % Gali (2008)

rrho_v       = 0.0;
rrho_a       = 0.9;
rrho_u       = 0.9;


ssigma_v = 0.1/100;
ssigma_a = 0.05/100;
ssigma_u = 0.05/100;

 



cal_para(1,1)  = bbeta;
cal_para(2,1)  = inv_ssigma;
cal_para(3,1)  = rrho_i;
cal_para(4,1)  = pphi_pi;
cal_para(5,1)  = pphi_ytilde;
cal_para(6,1)  = ppsi_yan;
cal_para(7,1)  = rrho_v;
cal_para(8,1)  = rrho_a;
cal_para(9,1) = rrho_u;
cal_para(10,1) = ssigma_v;
cal_para(11,1) = ssigma_a;
cal_para(12,1) = ssigma_u;


end