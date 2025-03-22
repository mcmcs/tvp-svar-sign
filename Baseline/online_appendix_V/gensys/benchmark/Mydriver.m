%%%%%%%%%%%%%%%%%
% 0. Housekeeping
%%%%%%%%%%%%%%%%%
clear variables; clc; close all; % clear global environment

% Declare auxiliary paths
addpath('helpfunctions')
addpath('helpfunctions/ChrisSimsGensys')




%%%%%%%%%%%%%%%%%%%%
% 1. Preliminaries %
%%%%%%%%%%%%%%%%%%%%

robust_dynare       = 1;       % Flag to check that Gensys policy functions equal Dynare policy functions
tol_dynare = 1e-5;

currdir=pwd;
%%%%%%%%%%%%%%%%%%%%%
% 2. Load our model %
%%%%%%%%%%%%%%%%%%%%%
% We linearize our model in  Mathematica using the file linearize_model.nb.
% The lineared model can be found is then obtained with the funciton
% get_model().


PC = 'steep_PC';


    switch PC
        case 'flat_PC'
        ttheta0 = 0.1;
        case 'steep_PC'
        ttheta0 = 0.5;
    end


[GAM0, GAM1, PPSI,PPI,C]            = get_model(ttheta0);           % Linearized model evaluated at ttheta0
[T1,C,T0,fmat,fwt,ywt,gev,eu,loose] = gensys(GAM0,GAM1,C,PPSI,PPI); % Policy functions evaluated at ttheta0

   
keyboard

if robust_dynare==1
 
    switch PC
        case 'flat_PC'
        load('dynare_flat_PC.mat');
        [out,errors,model_ss,dynare_ss] = get_dynare_check(tol_dynare,ttheta0,oo_flat);
        case 'steep_PC'
        load('dynare_steep_PC.mat');
        [out,errors,model_ss,dynare_ss] = get_dynare_check(tol_dynare,ttheta0,oo_steep);
    end
    
    if out==1
        display(['Steady State and Policy Functions equal to Dynare Solutions with a tolerance equal to ',num2str(tol_dynare)])
    else
        display(['Steady State and Policy Functions do not equal to Dynare Solutions: Check that the model is evaluated at the same parameters values'])
    end
    

end




% compute IRFs and compare with dynare
% y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.
nperiods = 21;
nshocks  = size(T0,2);
nvar = size(T0,1);
irf      = nan(nvar,nshocks,nperiods);

for h=1:nperiods


    irf(:,:,h) = (T1^(h-1))*T0;

end

 switch PC
        case 'flat_PC'
            gensys_flat.T1 = T1;
            gensys_flat.C = C;
            gensys_flat.T0 = T0;
            gensys_flat.fmat= fmat;
            gensys_flat.fwt=fwt;
            gensys_flat.ywt = ywt;
            gensys_flat.gev = gev;
            gensys_flat.eu = eu;
            gensys_flat.loose = loose;
            gensys_flat.irf = irf;
            save('gensys_flat_PC.mat','gensys_flat');

        case 'steep_PC'
            gensys_steep.T1 = T1;
            gensys_steep.C = C;
            gensys_steep.T0 = T0;
            gensys_steep.fmat= fmat;
            gensys_steep.fwt=fwt;
            gensys_steep.ywt = ywt;
            gensys_steep.gev = gev;
            gensys_steep.eu = eu;
            gensys_steep.loose = loose;
            gensys_steep.irf = irf;
            save('gensys_steep_PC.mat','gensys_steep');

 end








