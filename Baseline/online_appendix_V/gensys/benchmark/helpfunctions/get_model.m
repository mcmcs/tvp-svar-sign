function [GAM0,GAM1,PPSI,PPI,C] = get_model(ttheta)





calibrated_parameters = get_calibrated_parameters();

%% calibrated parameters
bbeta = calibrated_parameters(1,1);
inv_ssigma =calibrated_parameters(2,1);
rrho_i = calibrated_parameters(3,1);
pphi_pi = calibrated_parameters(4,1);
pphi_ytilde = calibrated_parameters(5,1);
ppsi_yan = calibrated_parameters(6,1);
rrho_v = calibrated_parameters(7,1);
rrho_a= calibrated_parameters(8,1);
rrho_u = calibrated_parameters(9,1);
ssigma_v = calibrated_parameters(10,1);
ssigma_a = calibrated_parameters(11,1);
ssigma_u =calibrated_parameters(12,1);


[model_ss,model_para_ss] = get_model_ss(ttheta);


kkappa      = ttheta(1);

ytilde = model_ss(1); 
ppi    = model_ss(2);
R      = model_ss(3);         
rn     = model_ss(4);
v      = model_ss(5);
a      = model_ss(6);
u      = model_ss(7);

rrho = model_para_ss(1);

% STATE VECTOR
% (1) = ytilde
% (2) = ppi
% (3) = R
% (4) = rn
% (5) = v
% (6) = a
% (7) = u
% (8) = E(ytilde(+1)
% (9)= E(ppi(t+1))
%varexo 

%vareps_v vareps_a vareps_u

nstate = 9;
nshocks=3;
nexp   = 2;
GAM0 =zeros(nstate,nstate);
    
GAM1 = zeros(nstate,nstate);
PPSI = zeros(nstate,nshocks);

PPI=zeros(nstate,nexp);
C=zeros(nstate,1);

%%//1. Euler equation
%ytilde = -inv_ssigma*(R-ppi(+1)-rn) + ytilde(+1);
GAM0(1,1) = 1;
GAM0(1,3) = inv_ssigma;
GAM0(1,9) = -inv_ssigma;
GAM0(1,4) = -inv_ssigma;
GAM0(1,8) = -1;


%//2. PC
% ppi = bbeta*ppi(+1) + kkappa*ytilde + u;
GAM0(2,2) = 1;
GAM0(2,9) = -bbeta;
GAM0(2,1) = -kkappa;
GAM0(2,7) = -1;


%//3. MP equation
% R = rrho_i*R(-1) + (1-rrho_i)*(rrho + pphi_pi*ppi +pphi_ytilde*ytilde) + v;
GAM0(3,3) = 1;
GAM0(3,2) = -(1-rrho_i)*pphi_pi;
GAM0(3,1) = -(1-rrho_i)*pphi_ytilde;
GAM0(3,5) = -1;
GAM1(3,3) = rrho_i;
C(3,1)    = (1-rrho_i)*rrho;

%//4. Natural rate of interest
%rn = rrho + (1/inv_ssigma)*ppsi_yan*(rrho_a-1)*a;
GAM0(4,4) = 1;
GAM0(4,6) = -(1/inv_ssigma)*ppsi_yan*(rrho_a-1);
C(4,1)    = rrho;


%//5. Monetary policy shock
%v=rrho_v*v(-1)  + ssigma_v*vareps_v;
GAM0(5,5) = 1;
GAM1(5,5) = rrho_v;

PPSI(5,1) = ssigma_v;

%//6. Technology shock
%a=rrho_a*a(-1)  + ssigma_a*vareps_a;
GAM0(6,6) = 1;
GAM1(6,6) = rrho_a;
PPSI(6,2) = ssigma_a;

%//7. Cost push shock
%u=rrho_u*u(-1)  + ssigma_u*vareps_u;
GAM0(7,7) = 1;
GAM1(7,7) = rrho_u;
PPSI(7,3) = ssigma_u;


% expectation error:
% ytilde (t) = E(t-1) ytilde(t) + \eta_y
GAM0(8,1) = 1;
GAM1(8,8) = 1;
PPI(8,1)  = 1;

% expectation error:
% ytilde (t) = E(t-1) ytilde(t) + \eta_y
GAM0(9,2) = 1;
GAM1(9,9) = 1;
PPI(9,2)  = 1;