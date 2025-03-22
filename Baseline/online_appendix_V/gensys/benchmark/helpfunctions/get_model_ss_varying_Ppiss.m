function [model_ss,model_para_ss] = get_model_ss_varying_Ppiss(ttheta,Ppiss)


calibrated_parameters = get_calibrated_parameters();

%% calibrated parameters
ssigmaa    = calibrated_parameters(1,1);  % Smets and Wouters (AER 2007)
ddelta     = calibrated_parameters(2,1);  % Steady state depreciation rate. Standard value in the literature
aalpha     = calibrated_parameters(3,1);  % Schmidt-Grohe and Uribe (ECTA 2012). This value is between Smets and Wouters (AER 2007) and Justiniano and Primiceri (2008)
eetaG      = calibrated_parameters(4,1);  % Average for our sample (Government Consumption Expenditures + Government Investment Expenditures)/GDP
iiotap     = calibrated_parameters(5,1);  % Christiano, Eichenbaum, and Trabandt (ECTA 2016)
iiotaw     = calibrated_parameters(6,1);  % Christiano, Eichenbaum, and Trabandt (ECTA 2016)
iiotaz     = calibrated_parameters(7,1);  % Christiano, Eichenbaum, and Trabandt (ECTA 2016)
%Ppiss      = calibrated_parameters(8,1);  % Implies an average annualized inflation rate equal to 2.46 percent (i.e. average GDP Deflator inflation for our sample)
Ztildess   = calibrated_parameters(9,1);  % Average real GDP growth in our sample 
mmutildess = calibrated_parameters(10,1); % Pc/Pi Average of our sample. Data obtained from Justiniano Primiceri and Tambalotti (RED 2011)


%% estimated parameters
% Monetary policy equation
% para 1  ppsi_ppi 
% para 2  ppsi_y  
% para 3  ppsi_gy
% para 4  rrho_R1
% para 5  rrho_R2

% Structural parameters
% para 6  kkappa
% para 7  b
% para 8  eeta_p
% para 9  eeta_w
% para 10 ttau
% para 11 bbeta
% para 12 nnu_p
% para 13 nnu_w
b           = ttheta(7);
eeta        = (1+ttheta(8))/ttheta(8);
eetaw       = (1+ttheta(9))/ttheta(9);
ttau        = ttheta(10);
bbeta       = 100/(100+ttheta(11)); 
nnu         = ttheta(12);
nnuw        = ttheta(13);

% Exogenous process: persistance
% para 14 rrho_d
% para 15 rrho_dL
% para 16 rrho_A
% para 17 rrho_mI
% para 18 rrho_mmu
% para 19 rrho_G


%% steady state and matrices for GENSYS
mI_ss=1;
Qtilde_ss = 1;
ggamma1 = Qtilde_ss*((Ztildess*mmutildess/bbeta)-1+ddelta);
ggamma2 = ggamma1*ssigmaa;
rtilde_ss = ggamma1;
rtildenopi_ss = ggamma1;
R_ss = Ztildess*Ppiss/bbeta;
Ppinopi_ss = 1;
Rnopi_ss = Ztildess*Ppinopi_ss/bbeta;
Ppistartilde_ss = ((1/(1-nnu))*(1-nnu*Ppiss^(-(1-iiotap)*(1-eeta))))^(1/(1-eeta));
Ppistartildenopi_ss = ((1/(1-nnu))*(1-nnu*Ppinopi_ss^(-(1-iiotap)*(1-eeta))))^(1/(1-eeta));
MCtilde1_ss = (1-nnu*bbeta*Ppiss^((1-iiotap)*eeta));
MCtilde1nopi_ss = (1-nnu*bbeta*Ppinopi_ss^((1-iiotap)*eeta));
MCtilde2_ss = (1-nnu*bbeta*Ppiss^(-(1-iiotap)*(1-eeta)));
MCtilde2nopi_ss = (1-nnu*bbeta*Ppinopi_ss^(-(1-iiotap)*(1-eeta)));
MCtilde_ss = ((eeta-1)/eeta)*(MCtilde1_ss/MCtilde2_ss)*Ppistartilde_ss;
MCtildenopi_ss = ((eeta-1)/eeta)*(MCtilde1nopi_ss/MCtilde2nopi_ss)*Ppistartildenopi_ss;
Wtilde_ss = (1/R_ss)*(1-aalpha)*(MCtilde_ss*(aalpha/rtilde_ss)^aalpha)^(1/(1-aalpha));
Wtildenopi_ss = (1/Rnopi_ss)*(1-aalpha)*(MCtildenopi_ss*(aalpha/rtildenopi_ss)^aalpha)^(1/(1-aalpha));
ratioWtildestarWtilde1_ss =(1-nnuw*(Ztildess^(-(1-eetaw)))*(Ppiss^(-(1-eetaw)*(1-iiotaw)))*(Ztildess^(iiotaz*(1-eetaw))));
ratioWtildestarWtilde1nopi_ss =(1-nnuw*(Ztildess^(-(1-eetaw)))*(Ppinopi_ss^(-(1-eetaw)*(1-iiotaw)))*(Ztildess^(iiotaz*(1-eetaw))));
ratioWtildestarWtilde2_ss = (1-nnuw);
ratioWtildestarWtilde2nopi_ss = (1-nnuw);
ratioWtildestarWtilde_ss = ((ratioWtildestarWtilde1_ss/ratioWtildestarWtilde2_ss)^(1/(1-eetaw)));
ratioWtildestarWtildenopi_ss = ((ratioWtildestarWtilde1nopi_ss/ratioWtildestarWtilde2nopi_ss)^(1/(1-eetaw)));
Wtildestartilde_ss = ratioWtildestarWtilde_ss*Wtilde_ss;Wstartilde_ss = Wtildestartilde_ss;
Wtildestartildenopi_ss = ratioWtildestarWtildenopi_ss*Wtildenopi_ss;
vp_ss = ((1-nnu)*Ppistartilde_ss^-eeta)/(1-nnu*Ppiss^((1-iiotap)*eeta));
vw1_ss = (1-nnuw);
vw2_ss = 1-nnuw*(Ppiss^((1-iiotaw)*eetaw))*(Ztildess^(iiotaz*(-eetaw)))*(Ztildess^eetaw);
vw_ss = (vw1_ss/vw2_ss)*(ratioWtildestarWtilde_ss^(-eetaw));
Oomega_ss = R_ss*(aalpha/(1-aalpha))*Wtilde_ss*Ztildess*mmutildess/rtilde_ss;
Oomeganopi_ss = Rnopi_ss*(aalpha/(1-aalpha))*Wtildenopi_ss*Ztildess*mmutildess/rtildenopi_ss;
aux_ss = ((1-eetaG)*MCtilde_ss*(Ztildess^-aalpha)*(mmutildess^-aalpha)*Oomega_ss^aalpha)-((Ztildess*mmutildess-1+ddelta)/(Ztildess*mmutildess))*Oomega_ss;
auxnopi_ss = ((1-eetaG)*MCtildenopi_ss*(Ztildess^-aalpha)*(mmutildess^-aalpha)*Oomeganopi_ss^aalpha)-((Ztildess*mmutildess-1+ddelta)/(Ztildess*mmutildess))*Oomeganopi_ss;
f1_ss = (1-bbeta*nnuw*(Ppiss^(-(1-iiotaw)*(1-eetaw)))*(Ztildess^(iiotaz*(1-eetaw)))*(Ztildess^(eetaw-1)))^-1;
f1nopi_ss = (1-bbeta*nnuw*(Ppinopi_ss^(-(1-iiotaw)*(1-eetaw)))*(Ztildess^(iiotaz*(1-eetaw)))*(Ztildess^(eetaw-1)))^-1;
f2_ss = ((eetaw-1)/eetaw)*(1/aux_ss)*(Wtildestartilde_ss^(1-eetaw))*Wtilde_ss^eetaw;
f2nopi_ss = ((eetaw-1)/eetaw)*(1/auxnopi_ss)*(Wtildestartildenopi_ss^(1-eetaw))*Wtildenopi_ss^eetaw;
f_ss = f1_ss*f2_ss*(Ztildess-bbeta*b)/(Ztildess-b);
fnopi_ss = f1nopi_ss*f2nopi_ss*(Ztildess-bbeta*b)/(Ztildess-b);
Ld1_ss = 1-bbeta*nnuw*((((Ppiss^iiotaw)*(Ztildess^iiotaz))/Ppiss)^(-eetaw*(1+ttau)))*Ztildess^(eetaw*(1+ttau));
Ld1nopi_ss = 1-bbeta*nnuw*((((Ppinopi_ss^iiotaw)*(Ztildess^iiotaz))/Ppinopi_ss)^(-eetaw*(1+ttau)))*Ztildess^(eetaw*(1+ttau));
Ld2_ss = (ratioWtildestarWtilde_ss)^(-eetaw*(1+ttau));
Ld2nopi_ss = (ratioWtildestarWtildenopi_ss)^(-eetaw*(1+ttau));
ppsiL = Ld1nopi_ss*fnopi_ss/Ld2nopi_ss;
Ld_ss = (Ld1_ss*f_ss/(ppsiL*Ld2_ss))^(1/(1+ttau));
Ktilde_ss = Oomega_ss*Ld_ss;
L_ss = vw_ss*Ld_ss;
Xtilde_ss = ((Ztildess*mmutildess-1+ddelta)/(Ztildess*mmutildess))*Ktilde_ss;
Ctilde_ss = aux_ss*Ld_ss;
Xxitilde_ss = (1/Ctilde_ss)*(Ztildess-bbeta*b)/(Ztildess-b);
YdGtilde_ss = (Ctilde_ss + Xtilde_ss)/(1-eetaG);
Ydtilde_ss = YdGtilde_ss;
U_ss=1;
g1_ss = (Xxitilde_ss*MCtilde_ss*Ydtilde_ss)*(1-nnu*bbeta*Ppiss^(eeta*(1-iiotap)))^-1;
g2_ss = (eeta/(eeta-1))*g1_ss;
d_ss = 1;dL_ss = 1;
Atilde_ss = (Ztildess^(1-aalpha))*(mmutildess)^(-aalpha);
Gtilde_ss = eetaG*YdGtilde_ss;
%flexible price
FPnnu = 0;
FPnnuw=0;
FPQtilde_ss = 1;
FPrtilde_ss = ggamma1;
FPrtildenopi_ss = ggamma1;
FPPpi_ss = Ppiss;
FPPpinopi_ss = 1;
FPR_ss = Ztildess*FPPpi_ss/bbeta;
FPRnopi_ss = Ztildess*FPPpinopi_ss/bbeta;
FPPpistartilde_ss = ((1/(1-FPnnu))*(1-FPnnu*FPPpi_ss^(-(1-iiotap)*(1-eeta))))^(1/(1-eeta));
FPPpistartildenopi_ss = ((1/(1-FPnnu))*(1-FPnnu*FPPpinopi_ss^(-(1-iiotap)*(1-eeta))))^(1/(1-eeta));
FPMCtilde1_ss = (1-FPnnu*bbeta*FPPpi_ss^((1-iiotap)*eeta));
FPMCtilde1nopi_ss = (1-FPnnu*bbeta*FPPpinopi_ss^((1-iiotap)*eeta));
FPMCtilde2_ss = (1-FPnnu*bbeta*FPPpi_ss^(-(1-iiotap)*(1-eeta)));
FPMCtilde2nopi_ss = (1-FPnnu*bbeta*FPPpinopi_ss^(-(1-iiotap)*(1-eeta)));
FPMCtilde_ss = ((eeta-1)/eeta)*(FPMCtilde1_ss/FPMCtilde2_ss)*FPPpistartilde_ss;
FPMCtildenopi_ss = ((eeta-1)/eeta)*(FPMCtilde1nopi_ss/FPMCtilde2nopi_ss)*FPPpistartildenopi_ss;
FPWtilde_ss = (1/FPR_ss)*(1-aalpha)*(FPMCtilde_ss*(aalpha/FPrtilde_ss)^aalpha)^(1/(1-aalpha));
FPWtildenopi_ss = (1/FPRnopi_ss)*(1-aalpha)*(FPMCtildenopi_ss*(aalpha/FPrtildenopi_ss)^aalpha)^(1/(1-aalpha));
FPratioWtildestarWtilde1_ss =(1-FPnnuw*(Ztildess^(-(1-eetaw)))*(FPPpi_ss^(-(1-eetaw)*(1-iiotaw)))*(Ztildess^(iiotaz*(1-eetaw))));
FPratioWtildestarWtilde1nopi_ss =(1-FPnnuw*(Ztildess^(-(1-eetaw)))*(FPPpinopi_ss^(-(1-eetaw)*(1-iiotaw)))*(Ztildess^(iiotaz*(1-eetaw))));
FPratioWtildestarWtilde2_ss = (1-FPnnuw);
FPratioWtildestarWtilde2nopi_ss = (1-FPnnuw);
FPratioWtildestarWtilde_ss = ((FPratioWtildestarWtilde1_ss/FPratioWtildestarWtilde2_ss)^(1/(1-eetaw)));
FPratioWtildestarWtildenopi_ss = ((FPratioWtildestarWtilde1nopi_ss/FPratioWtildestarWtilde2nopi_ss)^(1/(1-eetaw)));
FPWtildestartilde_ss = FPratioWtildestarWtilde_ss*FPWtilde_ss;
FPWtildestartildenopi_ss = FPratioWtildestarWtildenopi_ss*FPWtildenopi_ss;
FPOomega_ss = FPR_ss*(aalpha/(1-aalpha))*FPWtilde_ss*Ztildess*mmutildess/FPrtilde_ss;
FPOomeganopi_ss = FPRnopi_ss*(aalpha/(1-aalpha))*FPWtildenopi_ss*Ztildess*mmutildess/FPrtildenopi_ss;
FPaux_ss = ((1-eetaG)*FPMCtilde_ss*(Ztildess^-aalpha)*(mmutildess^-aalpha)*FPOomega_ss^aalpha)-((Ztildess*mmutildess-1+ddelta)/(Ztildess*mmutildess))*FPOomega_ss;
FPauxnopi_ss = ((1-eetaG)*FPMCtildenopi_ss*(Ztildess^-aalpha)*(mmutildess^-aalpha)*FPOomeganopi_ss^aalpha)-((Ztildess*mmutildess-1+ddelta)/(Ztildess*mmutildess))*FPOomeganopi_ss;
FPf1_ss = (1-bbeta*FPnnuw*(FPPpi_ss^(-(1-iiotaw)*(1-eetaw)))*(Ztildess^(iiotaz*(1-eetaw)))*(Ztildess^(eetaw-1)))^-1;
FPf1nopi_ss = (1-bbeta*FPnnuw*(FPPpinopi_ss^(-(1-iiotaw)*(1-eetaw)))*(Ztildess^(iiotaz*(1-eetaw)))*(Ztildess^(eetaw-1)))^-1;
FPf2_ss = ((eetaw-1)/eetaw)*(1/FPaux_ss)*(FPWtildestartilde_ss^(1-eetaw))*FPWtilde_ss^eetaw;
FPf2nopi_ss = ((eetaw-1)/eetaw)*(1/FPauxnopi_ss)*(FPWtildestartildenopi_ss^(1-eetaw))*FPWtildenopi_ss^eetaw;
FPf_ss = FPf1_ss*FPf2_ss*(Ztildess-bbeta*b)/(Ztildess-b);
FPfnopi_ss = FPf1nopi_ss*FPf2nopi_ss*(Ztildess-bbeta*b)/(Ztildess-b);
FPLd1_ss = 1-bbeta*FPnnuw*((((FPPpi_ss^iiotaw)*(Ztildess^iiotaz))/FPPpi_ss)^(-eetaw*(1+ttau)))*Ztildess^(eetaw*(1+ttau));
FPLd1nopi_ss = 1-bbeta*FPnnuw*((((FPPpinopi_ss^iiotaw)*(Ztildess^iiotaz))/FPPpinopi_ss)^(-eetaw*(1+ttau)))*Ztildess^(eetaw*(1+ttau));
FPLd2_ss = (FPratioWtildestarWtilde_ss)^(-eetaw*(1+ttau));
FPLd2nopi_ss = (FPratioWtildestarWtildenopi_ss)^(-eetaw*(1+ttau));
FPppsiL = FPLd1nopi_ss*FPfnopi_ss/FPLd2nopi_ss;
FPLd_ss = (FPLd1_ss*FPf_ss/(FPppsiL*FPLd2_ss))^(1/(1+ttau));
FPKtilde_ss = FPOomega_ss*FPLd_ss;
FPXtilde_ss = ((Ztildess*mmutildess-1+ddelta)/(Ztildess*mmutildess))*FPKtilde_ss;
FPCtilde_ss = FPaux_ss*FPLd_ss;
FPXxitilde_ss = (1/FPCtilde_ss)*(Ztildess-bbeta*b)/(Ztildess-b);
FPYdGtilde_ss = (FPCtilde_ss + FPXtilde_ss)/(1-eetaG);
FPYdtilde_ss = FPYdGtilde_ss;
FPGtilde_ss = eetaG*FPYdtilde_ss;
FPU_ss=1;


Wtilde      = log(Wtilde_ss);       % 1
model_ss(1) = Wtilde;       
Ctilde      = log(Ctilde_ss);       % 2
model_ss(2) = Ctilde;  
Ld          = log(Ld_ss);           % 3
model_ss(3) = Ld;  
L           = log(L_ss);            % 4
model_ss(4) = L;  
Qtilde      = log(Qtilde_ss);       % 5
model_ss(5) = Qtilde;  
Xtilde      = log(Xtilde_ss);       % 6 
model_ss(6) = Xtilde;  
Ztilde      = log(Ztildess);        % 7
model_ss(7) = Ztilde;  
rtilde      = log(rtilde_ss);       % 8
model_ss(8) = rtilde;  
U           = log(U_ss);            % 9
model_ss(9) = U;
Xxitilde    = log(Xxitilde_ss);     % 10
model_ss(10)= Xxitilde;
R           = log(R_ss);            % 11
model_ss(11)= R;
f           = log(f_ss);            % 12
model_ss(12)= f;
g1          = log(g1_ss);           % 13
model_ss(13)= g1;
g2          = log(g2_ss);           % 14
model_ss(14)= g2;
MCtilde     = log(MCtilde_ss);      % 15
model_ss(15)= MCtilde ;
Ydtilde     = log(Ydtilde_ss);      % 16
model_ss(16)= Ydtilde;
Wstartilde  = log(Wstartilde_ss);   % 17
model_ss(17)= Wstartilde;
Ppistar     = log(Ppistartilde_ss); % 18 
model_ss(18)= Ppistar;
Ktilde      = log(Ktilde_ss);       % 19
model_ss(19)= Ktilde;
Ppi         = log(Ppiss);           % 20
model_ss(20)= Ppi ;
vw          = log(vw_ss);           % 21
model_ss(21)= vw;
vp          = log(vp_ss);           % 22
model_ss(22)= vp;
d           = log(d_ss);            % 23
model_ss(23)= d;
dL          = log(dL_ss);           % 24
model_ss(24)= dL;
YdGtilde    = log(YdGtilde_ss);     % 25
model_ss(25)= YdGtilde;
Gtilde      = log(Gtilde_ss);       % 26
model_ss(26)= Gtilde;
Atilde      = log(Atilde_ss);       % 27
model_ss(27)= Atilde;
mmutilde    = log(mmutildess);      % 28
model_ss(28)= mmutilde;
mI          = log(mI_ss);           % 29
model_ss(29)= mI;
lagYdGtilde = YdGtilde;             % 30
model_ss(30)= lagYdGtilde;
lagR        = R;                    % 31
model_ss(31)= lagR;
lagCtilde   = Ctilde;               % 32
model_ss(32)= lagCtilde;
lagXtilde   = Xtilde;               % 33
model_ss(33)= lagXtilde;
lagWtilde   = Wtilde;               % 34
model_ss(34)= lagWtilde;
FPWtilde    = log(FPWtilde_ss);     % 35
model_ss(35)= FPWtilde;
FPCtilde    = log(FPCtilde_ss);     % 36
model_ss(36)= FPCtilde;
FPLd        = log(FPLd_ss);         % 37
model_ss(37)= FPLd;
FPQtilde    = log(FPQtilde_ss);     % 38
model_ss(38)= FPQtilde;
FPXtilde    = log(FPXtilde_ss);     % 39
model_ss(39)= FPXtilde;
FPrtilde    = log(FPrtilde_ss);     % 40
model_ss(40)= FPrtilde;
FPU         = log(FPU_ss);          % 41
model_ss(41)= FPU;
FPXxitilde  = log(FPXxitilde_ss);   % 42
model_ss(42)= FPXxitilde;
FPR         = log(FPR_ss);          % 43
model_ss(43)= FPR;
FPYdtilde   = log(FPYdtilde_ss);    % 44
model_ss(44)= FPYdtilde;
FPKtilde    = log(FPKtilde_ss);     % 45
model_ss(45)= FPKtilde;
FPPpi       = log(FPPpi_ss);        % 46
model_ss(46)= FPPpi;
FPYdGtilde  = log(FPYdGtilde_ss);   % 47
model_ss(47)= FPYdGtilde;
FPGtilde    = log(FPGtilde_ss);     % 48
model_ss(48)= FPGtilde;
FPlagR      = lagR;                 % 49
model_ss(49)= FPlagR;

model_para_ss(1,1) = ggamma1;
model_para_ss(2,1) = ggamma2;
model_para_ss(3,1) = ppsiL;
model_para_ss(4,1) = FPppsiL;


end
    
