% Algorithm that iterates
% Generating Q conditinal on B, C, D
% Generating DCD conditional on Q, B

% _jointpr: based on the joint prior, new single move sampler (06/25, 2024)

clc; close all;
clearvars -except ind_chain;

ind_chain = 3024 ;% 2024;


%% Data
curr_dir = pwd;
workpath = pwd;
frequency = 'quarterly';
data   = readtable([frequency,'_data/csvfiles/dataset.csv']);
data   = data(find(strcmp(string(data.dates),'1959-01-01')==1):find(strcmp(string(data.dates),'2023-04-01')==1),:);

data  = [100*diff(log(data.GDPC1)) 100*diff(log(data.PCEPILFE)) data.FEDFUNDS(2:end,1)/4 100*diff(log(data.M2)) data.BAA10YM(2:end,1)/4];

% addpath(genpath('dcd_helpfunctions'));
addpath(genpath('tvsvar_helpfunctions'));
addpath(genpath('var_helpfunctions'));

rng('default')
% rng(0)
seednum = 1234 * ind_chain;
rng(seednum)

%% Tuning parameters and restrictions

% -------------------------------------------------------------------------
% Main tuning parameters
scale_initial_V  = 4; % as in Primiceri, we scale the variance of the initial state variance by 4
scale_stepsize_V = 1; % control the overal size of the stepsize (relative to our original prior), sqrt(0.25)

% -------------------------------------------------------------------------
% Step size tuning parameters
k_W = 0.1 * scale_stepsize_V;
k_V = 0.1 * scale_stepsize_V;
k_B = 0.01 * scale_stepsize_V;
% -------------------------------------------------------------------------
% Restrictions:
type_restriction = 115; % 0 unrestricted, 1 baseline, 2 tyalor

switch type_restriction

    case 0
        function_restrictions = @restriction_none;
        is_restrict_B = false; %true if restriction includes B

    case 1
        function_restrictions = @restriction_baseline;
        is_restrict_B = false; %true if restriction includes B

    case 124
        function_restrictions = @restriction_baseline_noSys;
        is_restrict_B = false; %true if restriction includes B

    case 123
        function_restrictions = @restriction_baseline_noM;
        is_restrict_B = false; %true if restriction includes B

    case 111
        function_restrictions = @restriction_baseline_taylor;
        is_restrict_B = true; %true if restriction includes B

    case 112
        function_restrictions = @restriction_baseline_taylor_lg;
        is_restrict_B = true; %true if restriction includes B

    case 113
        function_restrictions = @restriction_baseline_taylor_1s;
        is_restrict_B = true; %true if restriction includes B

    case 114
        function_restrictions = @restriction_baseline_taylor_1s_gt0;
        is_restrict_B = true; %true if restriction includes B  

    case 115
        function_restrictions = @restriction_baseline_taylor_1s_gt0_gt1;
        is_restrict_B = true; %true if restriction includes B  

    case 11
        function_restrictions = @restriction_baseline_orig;
        is_restrict_B = false; %true if restriction includes B

    case 12
        function_restrictions = @restriction_baseline_orig_taylor;
        is_restrict_B = false; %true if restriction includes B

    case 2
        function_restrictions = @restriction_taylor;
        is_restrict_B = true; %true if restriction includes B
end

% -------------------------------------------------------------------------
% DCD_split
is_DCD_split = false; %true if new algorithm where D(t) and C(t) are separately updated (experimental)

%%
nburn = 0;
info_nvar = size(data,2);
T0 = 40;
lags =2;
y = data;
x = [];
for i=1:lags
    x=[x lag(y,i)];
end
x=[ones(length(x),1) x];
x=x(lags+1:end,:);
y=y(lags+1:end,:);
[T,n]=size(y);
x0=[x(1:T0,:)];
x=[x(T0+1:end,:)];
y0=y(1:T0,:);
y=y(T0+1:end,:);
Y0=reshape(y0,n*T0,1);
Y=reshape(y,n*(T-T0),1);
Z0=kron(eye(n),x0);
Z=kron(eye(n),x);

% defining the relevant dimensions
T=T-T0;         % length of estimation sample
k=n+lags*n^2;   % # of VAR coefficients
info.nvar = info_nvar;
info.ng = info_nvar*(info_nvar-1)/2;


%% Auxiliary
info_T    = T;

info_k    = k;
info_m    = info_k/info_nvar;
info_nlag = (info_m-1)/info_nvar;

% =========================================================================
% SETTING UP THE PRIORS
% time varying coefficients in B

% --- initialB is based on training sample
priorB   = olsblock(y0,x0);
Bbar   = reshape(priorB.bhatols,k,1)';
resB   = priorB.resols;
VBbar  = kron(resB'*resB/length(resB),eye(k/n)/priorB.XX);

info.kB   = k_B; % --- step size parameters
info.T0B  = 105; %T0+15;
info.VBbar = VBbar;

bbet_hat   = Bbar; % --- initial value parameters
bbet_hat_V = scale_initial_V * VBbar;
bbet_hat_Vinv = inv(bbet_hat_V);

% --- for ddelta and ggamma

Dn                       = duplication_matrix(info_nvar);
Dnplus                   = (Dn'*Dn)\Dn';
g                        = @(x)vechSsigma2ddelta1ggamma1(x,info);
Dg                       = NumericalDerivative(g,vech(priorB.sigmahatols));

iden                     = eye(size(Dg,1));
Zddelta1                 = iden(:,1:info_nvar);
Zggamma1                 = iden(:,info_nvar +1:end);
varvechSsigmahat         = 2*Dnplus*kron(priorB.sigmahatols,priorB.sigmahatols)*Dnplus';
varddelta1ggamma1hat_all = Dg*varvechSsigmahat*Dg';

% --- step size parameters
info.V_ddelta_1  = diag(Zddelta1'*(varddelta1ggamma1hat_all/T0)*Zddelta1);
info.k_w         = k_W;
info.vBar_w      = 2;
info.ddeltaBar_w = (info.k_w^2)*info.vBar_w*info.V_ddelta_1;
info.V_ggamma_1  = diag(Zggamma1'*(varddelta1ggamma1hat_all/T0)*Zggamma1);
info.k_v         = k_V;
info.vBar_v      = 2;
info.ddeltaBar_v = (info.k_v^2)*info.vBar_v*info.V_ggamma_1;

% --- initial value parameters
Dhat             = diag(diag(priorB.sigmahatols)).^(1/2);
Chat             = (Dhat\priorB.sigmahatols)/Dhat;

ddelta_hat = 2*log(diag(Dhat)); %mean for the initial period (Dhat is standard deviation; ddelta_hat is variance)
ddelta_hat_V = scale_initial_V * info.V_ddelta_1; %variance for the initial period
ddelta_hat_Vinv = eye(info.nvar) / diag(ddelta_hat_V); %precision for the initial period

logmChat   = logm(Chat); %mean for the initial period
ggamma_hat = logmChat(logical(tril(ones(info.nvar),-1)));
ggamma_hat_V = scale_initial_V * info.V_ggamma_1; %variance for the initial period
ggamma_hat_Vinv = eye(info.ng) / diag(ggamma_hat_V); %precision for the initial period


%% Setting the initial value
if type_restriction == 0
    % Drawing the initial VB,W,V from the prior
    V0B = VBbar*T0B*kB^2;
    VB_old = iwishrnd(V0B,T0B);

    Vd_old=nan(info.nvar,1);
    for i_d=1:info.nvar
        Vd_old(i_d) = 1/gamrnd(info.vBar_w/2, inv(info.ddeltaBar_w(i_d)/2));
    end

    Vg_old=nan(info.ng,1);
    for i_g=1:info.ng
        Vg_old(i_g) = 1/gamrnd(info.vBar_v/2, inv(info.ddeltaBar_v(i_g)/2));
    end

    % drawing the initial beta from the prior
    B_old(1,:) = reshape(priorB.bhatols,k,1)';
    for t=2:T
        B_old(t,:) =  reshape(priorB.bhatols,k,1)';
    end

    Dhat0  = diag(diag(priorB.sigmahatols)).^(1/2);
    Chat0  = (Dhat0\priorB.sigmahatols)/Dhat0;

    Ssigma0 = Dhat0*Chat0*Dhat0;
    condSsigma0 = cond( Ssigma0);

    ddelta_hat0 = 2*log(diag(Dhat0)); %mean for the initial period

    logmChat0   = logm(Chat0); %mean for the initial period
    ggamma_hat0 = logmChat0(logical(tril(ones(info.nvar),-1)));

    % drawing the initial delta from the prior
    ddelta_old(1,:) =     ddelta_hat0;
    for t=2:T
        ddelta_old(t,:) =  ddelta_hat0;
    end


    % drawing the initial ggamma from the prior
    ggamma_old(1,:) = ggamma_hat0;
    for t=2:T
        ggamma_old(t,:) = ggamma_hat0;
    end

    DCDhalf_old = nan(info_nvar, info_nvar, info_T);
    for t=1:info_T
        D_old = diag(exp(ddelta_old(t,:)/2));
        C_old = ggammatoC(ggamma_old(t,:));
        Ssigma_old = D_old*C_old*D_old;
        DCDhalf_old(:,:,t) = chol(Ssigma_old, 'lower');
    end


else
    %% Posterior mean
    %load(['results',filesep,'temp_results_unrest_ndraws_50000_is_restrict_B0.mat'], 'mat_B', 'mat_ddelta', 'mat_ggamma', 'mat_Vd', 'mat_Vg', 'mat_VB');
     load(['results',filesep,'temp_results_jointpr_rest115r2024ndraws500000.mat'], 'mat_B', 'mat_ddelta', 'mat_ggamma', 'mat_Vd', 'mat_Vg', 'mat_VB');

    B_old = mat_B(:,:,end);
    ddelta_old  = mat_ddelta(:,:,end);
    ggamma_old  = mat_ggamma(:,:,end);
    DCDhalf_old = nan(info_nvar, info_nvar, info_T);
    for t=1:info_T
        D_old = diag(exp(ddelta_old(t,:)/2));
        C_old = ggammatoC(ggamma_old(t,:));
        Ssigma_old = D_old*C_old*D_old;
        DCDhalf_old(:,:,t) = chol(Ssigma_old, 'lower');
    end
    % Initial W_old and V_old
    Vd_old = mat_Vd(end,:);
    Vg_old = mat_Vg(end,:);
    VB_old = mat_VB(:,:,end);


    %     load('temp_results_unrest_ndraws_50000_is_restrict_B0.mat');
    %     B_old = mean(mat_B(:,:,nburn+1:end),3);
    %     ddelta_old  = mean(mat_ddelta(:,:,nburn+1:end),3);
    %     ggamma_old  = mean(mat_ggamma(:,:,nburn+1:end),3);
    %     DCDhalf_old = nan(info_nvar, info_nvar, info_T);
    %     for t=1:info_T
    %         D_old = diag(exp(ddelta_old(t,:)/2));
    %         C_old = ggammatoC(ggamma_old(t,:));
    %         Ssigma_old = D_old*C_old*D_old;
    %         DCDhalf_old(:,:,t) = chol(Ssigma_old, 'lower');
    %     end

    %     %% Initial W_old and V_old
    %     Vd_old = mean(mat_Vd(nburn+1:end,:),1);
    %     Vg_old = mean(mat_Vg(nburn+1:end,:),1);
    %     VB_old = mean(mat_VB(:,:,nburn+1:end),3);

    %% Initial Q
    Q_old = nan(info_nvar, info_nvar, info_T);

    for t=1:info_T
        Q_old(:,:,t) = draw_from_On(info_nvar);
    end

end

%% Iteration - multiple runs
% preliminaries
Yhat = (Y-sum((Z.*repmat(squeeze(B_old),info_nvar,1))')');
yhat = reshape(Yhat,info_T,info_nvar); % residual

% setting
M0           = 1000000;          % number of MCMC draws
save_every   = 50;             % save every save_every
report_every = 10;          % report every report_every
nsave      = M0/save_every; % we store nsave elements
ind_save   = 1;             % counter for save

% --- matrix to store

mat_L0 = zeros(info_nvar,info_nvar,T,nsave);
mat_B  = zeros(T,k,nsave);
mat_Q  = zeros(info_nvar,info_nvar,T,nsave);

mat_ddelta = zeros(T,info_nvar,nsave);
mat_ggamma = zeros(T,info_nvar*(info_nvar-1)/2,nsave);

mat_VB = zeros(k,k,nsave);
mat_Vd = zeros(nsave,info_nvar);
mat_Vg = zeros(nsave,info_nvar*(info_nvar-1)/2);
mat_logL = zeros(nsave,1);
logL = zeros(M0,1);

L0_old = zeros(info_nvar); %matrix to initialize

mat_St = zeros(M0,1);


% To make sure that we start with a draw that satisifies restrictions
tic
disp('check 1')
parfor t=1:T
    Q_old(:,:,t) = draw_Q_condi(B_old(t,:), ddelta_old(t,:), ggamma_old(t,:), function_restrictions, info_nvar, t);
end
toc
disp('check 2')


% actual run
mean_L0_11_210 = 0; % value to monitor IRFs
time_to_finish = tic; %start clock
for i = 1:M0

    % =====================================================================
    % Drawing B
    % =====================================================================

    if is_restrict_B
        % if restriction includes B
        [B_old, ntry_B_old] = draw_B_single_move(y, Z, B_old, ddelta_old, ggamma_old, Q_old, VB_old, bbet_hat, bbet_hat_Vinv, function_restrictions);

    else
        % if restriction does not include B, we can do a standard "multi-move" sampling
        % Kalman filter
        SHAT=zeros(T,k);
        SIG=zeros(T,k,k);
        %sig=4*VBbar;
        %shat=Bbar';
        sig = bbet_hat_V;
        shat = bbet_hat';
        for t=1:T
            Ssigma_half = DCDhalf_old(:,:,t);
            Ssigma = Ssigma_half * Ssigma_half';
            %[shat,sig] = kfilter(y(t,:)',Z([0:n-1]*T+t,:),eye(k),shat,sig,Ssigma,squeeze(VB_old(1:k,1:k)));
            [shat,sig] = kfilter_JoE(y(t,:)',Z([0:n-1]*T+t,:),eye(k),shat,sig,Ssigma,squeeze(VB_old(1:k,1:k)),t); % in this way, the first observation does not depend on the variance of the step size
            SHAT(t,:)=shat';
            SIG(t,:,:)=sig;
        end

        % simulation smoother
        B_old(T,:)=mvnrnd(shat,sig/2+sig'/2,1);
        for t=T-1:-1:1
            [btTp,StTp]=kback(SHAT(t,:),squeeze(SIG(t,:,:)),squeeze(B_old(t+1,:)),eye(k),squeeze(VB_old(1:k,1:k)));
            B_old(t,:)=mvnrnd(btTp,StTp/2+StTp'/2,1);
        end
    end

    % update residual
    Yhat = (Y-sum((Z.*repmat(squeeze(B_old),info_nvar,1))')');
    yhat = reshape(Yhat,info_T,info_nvar); % residual


    % =====================================================================
    % Drawing Q
    %     tic
    %disp('Draw Q: begin')
    parfor t=1:T
        Q_old(:,:,t) = draw_Q_condi(B_old(t,:), ddelta_old(t,:), ggamma_old(t,:), function_restrictions, info_nvar, t)
    end
    %disp('Draw Q: end')
    %     toc

    % =====================================================================
    % Drawing DCD
    %     tic;
    if ~is_DCD_split
        [ddelta_old, ggamma_old, ntry_old] = draw_DCD_single_move(yhat, B_old, ddelta_old, ggamma_old, Q_old, Vd_old, Vg_old, ddelta_hat, ddelta_hat_Vinv, ggamma_hat, ggamma_hat_Vinv, function_restrictions);
    else
        [ddelta_old, ggamma_old, ntry_old] = draw_DCD_single_move_split(yhat, B_old, ddelta_old, ggamma_old, Q_old, Vd_old, Vg_old, ddelta_hat, ddelta_hat_Vinv, ggamma_hat, ggamma_hat_Vinv, function_restrictions);
    end

    for t=1:info_T
        Dt_old    = diag(exp(ddelta_old(t,:)/2));
        Ct_old    = ggammatoC(ggamma_old(t,:));
        DCDhalf_old(:,:,t) =  chol(Dt_old*Ct_old*Dt_old, 'lower');
    end
    %     toc;

    % =====================================================================
    % Drawing var(errB)
    %     kB  = 0.01;
    %     T0B = 55;
    V0B = info.VBbar*info.T0B*info.kB^2;

    errB=[B_old(2:end,:)-B_old(1:end-1,:)];
    V1B=errB'*errB+V0B;
    VB_old=iwishrnd(V1B,T-1+info.T0B);

    % =====================================================================
    % STEP 5: DRAWS OF Vd
    Vd_prop = Vd_old;
    errD = ddelta_old(2:end,:) - ddelta_old(1:end-1,:);
    vTilde_w = size(errD,1) + info.vBar_w;
    for i_ig=1:info.nvar
        ddeltaTilde_w = info.ddeltaBar_w(i_ig) + sum(errD(:,i_ig).^2);
        Vd_prop(i_ig) = 1/gamrnd(vTilde_w/2, inv(ddeltaTilde_w/2));
    end
    Vd_old = Vd_prop; % we don't need MH correction for the first step because the distribution of the initial delta does not depend on Vd

    % =====================================================================
    % STEP 6: DRAWS OF Vg
    Vg_prop = Vg_old;
    errG = ggamma_old(2:end,:) - ggamma_old(1:end-1,:);
    vTilde_v = size(errG,1) + info.vBar_v;
    for i_ig=1:info.ng
        ddeltaTilde_v = info.ddeltaBar_v(i_ig) + sum(errG(:,i_ig).^2);
        Vg_prop(i_ig)  = 1/gamrnd(vTilde_v/2, inv(ddeltaTilde_v/2));
    end
    Vg_old = Vg_prop; % we don't need MH correction for the first step because the distribution of the initial gamma does not depend on Vg

    % =====================================================================
    % IRF on impact
    for t=1:info_T
        L0_old(:,:,t) = DCDhalf_old(:,:,t) * Q_old(:,:,t);
    end
    mean_L0_11_210 = mean_L0_11_210 * (i-1)/i + L0_old(1,1,210)/i;

    % =====================================================================
    % likelihood calculation (conditional on B and Sigma)
    temp_lik = zeros(T,1);
    for t=1:T
        temp_lik(t,:) = log( mvnpdf(y(t,:)', Z([0:n-1]*T+t,:) * B_old(t,:)', DCDhalf_old(:,:,t) * DCDhalf_old(:,:,t)') );
    end
    logL(i,:) = sum(temp_lik); %log-like of i's draw


    % =====================================================================
    % report/store
    if rem(i, report_every) == 0
        % --- report
        disp(' ====================================== ');
        disp([num2str(i), '-th draw / ', num2str(M0)]);
        ttt = toc(time_to_finish);
        disp(['remained time (sec.) = ', num2str( (ttt/i) * (M0-i) ) ]);
        disp(['time per draw (sec.) = ', num2str( (ttt/i) ) ]);
        disp(['loglik = ', num2str(logL(i,1))]);
        disp(['L0 = ', num2str(L0_old(1,1,210))]);
        disp(['mean L0 = ', num2str( mean_L0_11_210 )]);
    end

    if rem(i, save_every) == 0
        % --- store
        mat_L0(:,:,:,ind_save)   = L0_old;
        mat_B(:,:,ind_save)      = B_old;
        mat_Q(:,:,:,ind_save)    = Q_old;
        mat_ddelta(:,:,ind_save) = ddelta_old;
        mat_ggamma(:,:,ind_save) = ggamma_old;
        mat_VB(:,:,ind_save)     = VB_old;
        mat_Vd(ind_save,:)       = Vd_old;
        mat_Vg(ind_save,:)       = Vg_old;
        mat_logL(ind_save,1)     = logL(i,1);

        ind_save = ind_save + 1;
    end

end

% Store temporary results
switch type_restriction

    case 0
        save(['temp_results_unrest','_ndraws_',num2str(M0),'_is_restrict_B',num2str(is_restrict_B),'.mat'],'mat_B','mat_ddelta','mat_ggamma','mat_VB','mat_Vd','mat_Vg','mat_logL');
    otherwise
        cd results
        % save(['temp_results_jointpr_rest', num2str(type_restriction), 'r' ,num2str(ind_chain),'ndraws_',num2str(M0),'is_restrict_B',num2str(is_restrict_B),'.mat'],'mat_logL','mat_B','mat_L0','y','x');
        save(['temp_results_jointpr_rest', num2str(type_restriction), 'r' ,num2str(ind_chain),'ndraws',num2str(M0),'.mat']);
        cd ..

end

