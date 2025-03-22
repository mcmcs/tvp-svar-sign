clear variables; clc; close all; % clear global environment

% Declare auxiliary paths
addpath('helpfunctions')
addpath('helpfunctions/ChrisSimsGensys')


% compute IRFs and compare with dynare
% y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.

load('gensys_flat_PC.mat','gensys_flat');
load('gensys_steep_PC.mat','gensys_steep');


%% simulation with structural break
T_sim      = 215+42;
T_break    = 201+42;
nshocks = size(gensys_flat.T0,2);

simdata     = nan(size(gensys_flat.T0,1),T_sim);

% initial values for simulation
  ttheta0 = 0.1;
[model_ss,~] = get_model_ss(ttheta0);

simdata(:,1)     = [model_ss';zeros(2,1)];% last 2 zeros are expectational erros

for tt=2:T_sim
    if tt<=T_break
        simdata(:,tt) = gensys_flat.T1*simdata(:,tt-1) + gensys_flat.C + gensys_flat.T0*randn(nshocks,1);
    else
        simdata(:,tt) = gensys_steep.T1*simdata(:,tt-1) + gensys_steep.C + gensys_steep.T0*randn(nshocks,1);
    end
end


ytildeOBS = simdata(1,:)*100;
ppiOBS    = simdata(2,:)*400 + 2;
ROBS      = simdata(3,:)*400;

save('simdata.mat','simdata');

close all
figure(1)
subplot(1,3,1)
plot(ytildeOBS)
title('Output Gap (%)')
hold on
subplot(1,3,2)
plot(ppiOBS)
title('Inflation (% Annualized)')
hold on
subplot(1,3,3)
plot(ROBS)
title('FFR (% Annualized)')
