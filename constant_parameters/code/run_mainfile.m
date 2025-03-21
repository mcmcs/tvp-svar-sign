%==========================================================================
%% housekeeping
%==========================================================================
clear variables;
close all;
userpath('clear');
clc;

rng('default'); % reinitialize the random number generator to its startup configuration
rng(0);         % set seed

currdir=pwd;
cd ..
get_help_dir_currdir=pwd;
addpath([get_help_dir_currdir,'/helpfunctions']); % set path to helper functions
cd(currdir)

%==========================================================================
%% load the data and priors
%==========================================================================
%% Data
curr_dir = pwd;
workpath = pwd;
frequency = 'quarterly';
data   = readtable([frequency,'_data/csvfiles/dataset.csv']);
data   = data(find(strcmp(string(data.dates),'1959-01-01')==1):find(strcmp(string(data.dates),'2023-04-01')==1),:);
%data   = data(find(strcmp(string(data.dates),'1959-01-01')==1):find(strcmp(string(data.dates),'2007-04-01')==1),:);
%% data
% Minchul thinks that FFR also has to be QoQ growth/
num  = [100*diff(log(data.GDPC1)) 100*diff(log(data.PCEPILFE)) data.FEDFUNDS(2:end,1)/4 100*diff(log(data.M2)) data.BAA10YM(2:end,1)/4];


%==========================================================================
%% model setup
%==========================================================================
nlag      = 4;                  % number of lags
nvar      = size(num,2);                   % number of endogenous variables
nex       = 1;                   % set equal to 1 if a constant is included; 0 otherwise
m         = nvar*nlag + nex;     % number of exogenous variables
nd        = 1e5;                 % number of orthogonal-reduced-form (B,Sigma,Q) draws
iter_show = 1e4;                 % display iteration every iter_show draws
horizon   = 20;                  % maximum horizon for IRFs
index     = [4 8 16 24 32 40];                  % define  horizons for the FEVD
horizons  = 0:0;                 % horizons to restrict
NS        = 1 + numel(horizons); % number of objects in F(A_{0},A_{+}) to which we impose sign and zero restrictios: F(THETA)=[A_{0};L_{0},...,L_{horizons}]
e         = eye(nvar);           % create identity matrix
maxdraws  = 1e4;                 % max number of importance sampling draws
conjugate = 'none';        % structural or irfs or empty

%==========================================================================
%% identification: declare Ss and Zs matrices
%==========================================================================
% restrictions on A0  and/or IRFs

% sign restrictions
S = cell(nvar,1);
for ii=1:nvar
    S{ii}=zeros(0,nvar*NS);
end
ns1  = 3;
S{1} = zeros(ns1,nvar*NS);
S{1}(1,nvar+2)   =  -1;
S{1}(2,nvar+3)   =   1;
S{1}(3,nvar+4)   =  -1;

% % zero restrictions
Z=cell(nvar,1);
for i=1:nvar
    Z{i}=zeros(0,nvar*NS);
end

% nz1  = 2;
% Z{1} = zeros(nz1,nvar*NS);
% Z{1}(1,4)=1;
% Z{1}(2,5)=1;


%==========================================================================
%% Setup info
%==========================================================================
info=SetupInfo(nvar,m,Z,@(x)chol(x));

% ZF(A_{0},A_{+})
info.nlag     = nlag;
info.horizons = horizons;
info.ZF       = @(x,y)ZF(x,y);

% functions useful to compute the importance sampler weights
iw_info = info;
fs      = @(x)ff_h(x,iw_info);
r       = @(x)ZeroRestrictions(x,iw_info);

if strcmp(conjugate,'irfs')==1
    fo              = @(x)f_h(x,iw_info);
    fo_str2irfs     = @(x)StructuralToIRF(x,iw_info);
    fo_str2irfs_inv = @(x)IRFToStructural(x,iw_info);
    r_irfs          = @(x)IRFRestrictions_more_general(x,iw_info); 
end


% function useful to check the sign restrictions
fh_S_restrictions  = @(x)SF(x,iw_info,S);

%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [yt(t-1), ... , yt(t-nlag), constant];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
yt = num(nlag+1:end,:);
T  = size(yt,1);
xt = zeros(T,nvar*nlag+nex);
for i=1:nlag
    xt(:,nvar*(i-1)+1:nvar*i) = num((nlag-(i-1)):end-i,:) ;
end
if nex>=1
    xt(:,nvar*nlag+nex)=ones(T,1);
end
% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors
keyboard

%% prior for reduced-form parameters
nnuBar              = 0;
OomegaBarInverse    = zeros(m);
PpsiBar             = zeros(m,nvar);
PphiBar             = zeros(nvar);



%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = OomegaTilde*(X'*Y + OomegaBarInverse*PpsiBar);
PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;




%% useful definitions
% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws         = cell([nd,1]); % reduced-form lag parameters
Sigmadraws     = cell([nd,1]); % reduced-form covariance matrices
Qdraws         = cell([nd,1]); % orthogonal matrices
storevefh      = zeros(nd,1);  % volume element f_{h}
storevegfhZ    = zeros(nd,1);  % volume element g o f_{h}|Z
uw             = zeros(nd,1);  % unnormalized importance sampler weights

if strcmp(conjugate,'irfs')==1
    storevephi      = zeros(nd,1);  % volume element f_{h}
    storevegphiZ    = zeros(nd,1);  % volume element g o f_{h}|Z
end

% definitions related to IRFs; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(nvar),nlag-1,1)];
A      = cell(nlag,1);
extraF = repmat(zeros(nvar),1,nlag-1);
F      = zeros(nlag*nvar,nlag*nvar);
for l=1:nlag-1
    F((l-1)*nvar+1:l*nvar,nvar+1:nlag*nvar)=[repmat(zeros(nvar),1,l-1) e repmat(zeros(nvar),1,nlag-(l+1))];
end

% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below


%% initialize counters to track the state of the computations
counter = 1;
record  = 1;
count   = 0;
tStart = tic;
while record<=nd
    
    
    %% step 1 in Algorithm 2
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(PpsiTilde,nvar*m,1);
    Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    % store reduced-form draws
    Bdraws{record,1}     = Bdraw;
    Sigmadraws{record,1} = Sigmadraw;
    
   
    %% steps 2:4 of Algorithm 2
    w           = DrawW(iw_info);   
    x           = [vec(Bdraw); vec(Sigmadraw); w];
    structpara  = ff_h_inv(x,iw_info);
    
    % store the matrix Q associated with step 3
    Qdraw            = SpheresToQ(w,iw_info,Bdraw,Sigmadraw);
    Qdraws{record,1} = reshape(Qdraw,nvar,nvar);

    A0 = reshape(structpara(1:nvar^2,1),nvar,nvar);

    ppsiy = -A0(1,1)/A0(3,1);
    ppsip = -A0(2,1)/A0(3,1);
    ppsim = -A0(4,1)/A0(3,1);
    ppsics = -A0(5,1)/A0(3,1);


% 
%     non_stationary = 0;
%         Fmatrix = zeros(nvar*nlag,nvar*nlag);
%         for ell=1:nlag
%             F(1:nvar,1 + (ell-1)*nvar:ell*nvar)=Bdraw(1+(ell-1)*nvar:ell*nvar,:)';
%             if ell<nlag
%                 Fmatrix(1+ell*nvar:nvar+ell*nvar,1+(ell-1)*nvar:ell*nvar)=eye(nvar);
%             end
%         end
%         EF=eig(Fmatrix);
%         for i=1:nvar*nlag
%             if norm(EF(i))>1
%                 non_stationary = 1;
%             end
%         end
%  
   
    %% check if sign restrictions hold
    signs      = fh_S_restrictions(structpara);
    
    
    if (sum(signs{1}*e(:,1)>0))==size(signs{1}*e(:,1),1) %&& ppsiy>0 && ppsip>0 && ppsim>0 && ppsics<0 && ppsiy<4 && ppsip<4 && ppsim<4 && ppsics>-4 
        
        count=count+1;
        
        %% compute importance sampling weights
       
        switch conjugate
            
            case 'structural'
                
                
                storevefh(record,1)   = (nvar*(nvar+1)/2)*log(2)-(2*nvar+m+1)*LogAbsDet(reshape(structpara(1:nvar*nvar),nvar,nvar));
                storevegfhZ(record,1) = LogVolumeElement(fs,structpara,r); 
                uw(record,1)          = exp(storevefh(record,1) - storevegfhZ(record,1));
                
            case 'irfs'
                
                irfpara                = fo_str2irfs(structpara);
                storevephi(record,1)   = LogVolumeElement(fo,structpara)   + LogVolumeElement(fo_str2irfs_inv,irfpara);
                storevegphiZ(record,1) = LogVolumeElement(fs,structpara,r) + LogVolumeElement(fo_str2irfs_inv,irfpara,r_irfs); 
                uw(record,1)           = exp(storevephi(record,1) - storevegphiZ(record,1));
                
            otherwise
                
                uw(record,1) = 1;
                
        end
        
    else
        
        uw(record,1) = 0;
        
    end
    
    if counter==iter_show
        
        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(nd-(record))])
        counter =0;
        
    end
    counter = counter + 1;
    record=record+1;
    
end


tElapsed = toc(tStart);
imp_w  = uw/sum(uw);
ne = floor(1/sum(imp_w.^2));


%% store draws
Ltilde        = zeros(horizon+1,nvar,nvar,ne); % define array to store IRF
CLtilde        = zeros(horizon+1,nvar,nvar,ne); % define array to store IRF
A0tilde       = zeros(nvar,nvar,ne); % define array to store A0
Aplustilde    = zeros(m,nvar,ne); % define array to store Aplus
FEVD          = zeros(nvar,size(index,2),ne);                 % define array to store FEVD
hist_is_draws = zeros(ne,1);   % define array to store draws from importance sampler



for s=1:min(ne,maxdraws)
    
    %% draw: B,Sigma,Q
    is_draw = randsample(1:size(imp_w,1),1,true,imp_w);
    hist_is_draws(s,1)=is_draw;

    Bdraw       = Bdraws{is_draw,1};
    Sigmadraw   = Sigmadraws{is_draw,1};
    Qdraw       = Qdraws{is_draw,1};
    
    
    % compute matrix F: useful for FEVD
    hSigmadraw = hh(Sigmadraw);
    A0         = hSigmadraw\e;
    Aplus      = Bdraw*A0;
    % Page 8 ARRW
    for l=1:nlag-1
        A{l} = Aplus((l-1)*nvar+1:l*nvar,1:end);
        F((l-1)*nvar+1:l*nvar,1:nvar)=A{l}/A0;
    end
    A{nlag} = Aplus((nlag-1)*nvar+1:nlag*nvar,1:end);
    F((nlag-1)*nvar+1:nlag*nvar,:)=[A{nlag}/A0 extraF];
    
    x=[reshape(Bdraw,m*nvar,1); reshape(Sigmadraw,nvar*nvar,1); Qdraw(:)];
    structpara = f_h_inv(x,info);
    
    
    LIRF =IRF_horizons(structpara, nvar, nlag, m, 0:horizon);
    
    
    for h=0:horizon
        Ltilde(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);

        for i_shock=1:info.nvar
                CLtilde(h+1,[1 2 4],i_shock,s)   = sum(Ltilde(1:h+1,[1 2 4],i_shock,s),1);
                CLtilde(h+1,[3 5],i_shock,s)       = Ltilde(h+1,[3 5],i_shock,s);
            end

    end
    A0tilde(:,:,s)    = reshape(structpara(1:nvar*nvar),nvar,nvar);
    Aplustilde(:,:,s) = reshape(structpara(nvar*nvar+1:end),m,nvar);
    
    
    
end

A0tilde    = A0tilde(:,:,1:s);
Aplustilde = Aplustilde(:,:,1:s);
Ltilde     = Ltilde(:,:,:,1:s);
CLtilde     = CLtilde(:,:,:,1:s);

cd results
savefile='results.mat';
save(savefile,'Ltilde','CLtilde','A0tilde','Aplustilde','imp_w','ne','Y','X','nlag');
cd ..



