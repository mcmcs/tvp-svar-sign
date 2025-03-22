%% Housekeeping
clc; clear variables; close all;


%% Load results 
current_path = pwd;
addpath(genpath('figures_helpfunctions'));
addpath(genpath('tvsvar_helpfunctions'));
addpath(genpath('var_helpfunctions'));
results_path = pwd;
cd(current_path);


specname = 'rest115r3024ndraws1000000';

load(['results', filesep, ['temp_results_jointpr_',specname, '.mat']]);


mkdir(['results', filesep, specname]);
savepath = [pwd, filesep, 'results', filesep, specname];


%% burn 1/4 of mcmc chain
nburn00 = 5000;
max_draws = size(mat_Q,4);
mat_Q = mat_Q(:,:,:,nburn00+1:max_draws);
mat_B = mat_B(:,:,nburn00+1:max_draws);
mat_ddelta = mat_ddelta(:,:,nburn00+1:max_draws);
mat_ggamma = mat_ggamma(:,:,nburn00+1:max_draws);
mat_L0 = mat_L0(:,:,:,nburn00+1:max_draws);



%% Computing objects of interest for baseline identification scheme
%  1.1. Get the structural parameters
%  1.2. Compute the systematic component of monetary policy
%  1.3. Compute the historical decomposition
ndraws  = size(mat_L0,4);
A_old   = zeros(info_nvar);
e       = eye(info_nvar);

r2.Ssigma    = zeros(info_nvar,info_nvar,ndraws);
info.hbar    = 20; % max IRF horizon
Ltilde_RC    = nan(info.hbar+1,info_nvar,info_nvar,ndraws,info_T);
CLtilde_RC   = nan(info.hbar+1,info_nvar,info_nvar,ndraws,info_T);
shocks_RC    = nan(info_nvar,ndraws,info_T);
mcal_A_rw    = nan(info_nvar*info_nvar,info_T,ndraws);
mcal_F_rw    = nan(info_k,info_T,ndraws);


ttime                        = 1969.75:0.25:2023.25;

% Historical decomposition starting in 2022Q1
t_origin_2022Q1               = 2022;
index_t_origin_2022Q1        = find(ttime==t_origin_2022Q1);
Hhistdecomp_2022Q1           = 4;
fcast_2022Q1                 = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
shock_cont_2022Q1            = nan(info_nvar,info_nvar,Hhistdecomp_2022Q1+1,ndraws);
shocks_RC_fixed_para_2022Q1  = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);

% Historical decomposition starting in 2021Q1
t_origin_2021Q1              = 2021;
index_t_origin_2021Q1        = find(ttime==t_origin_2021Q1);
Hhistdecomp_2021Q1           = 8;
fcast_2021Q1                 = nan(info_nvar,Hhistdecomp_2021Q1+1,ndraws);
shock_cont_2021Q1            = nan(info_nvar,info_nvar,Hhistdecomp_2021Q1+1,ndraws);
shocks_RC_fixed_para_2021Q1  = nan(info_nvar,Hhistdecomp_2021Q1+1,ndraws);


for i = 1:ndraws

    disp([num2str(i), ' / ', num2str(ndraws)]);
    
    % impulse responses and structural parameters
    for t=1:info_T
        A_old(:,:,t)  = inv(squeeze(mat_L0(:,:,t,i)))';
        A = squeeze(A_old(:,:,t));
        B = reshape(mat_B(t,:,i),info_m,info_nvar);
        F=B*A;
        structpara = [A(:);F(:)];
        shocks_RC(:,i,t) = ((y(t,:)-x(t,:)*B)*A)';
        LIRF       = IRF_horizons(structpara, info_nvar, info_nlag, info_m, 0:info.hbar);
        for h=0:info.hbar
            
            Ltilde_RC(h+1,:,:,i,t) =  LIRF(1+h*info_nvar:(h+1)*info_nvar,:);
            for i_shock=1:info_nvar
                CLtilde_RC(h+1,[1 2 4],i_shock,i,t)   = sum(Ltilde_RC(1:h+1,[1 2 4],i_shock,i,t),1);
                CLtilde_RC(h+1,[3 5],i_shock,i,t)     = Ltilde_RC(h+1,[3 5],i_shock,i,t);
            end
        end
        mcal_A_rw(:,t,i) = A(:);
        mcal_F_rw(:,t,i) = F(:);
    end



    % historical decomposition 2022Q1
    A          = reshape(mcal_A_rw(:,index_t_origin_2022Q1,i),info_nvar,info_nvar);
    F          = reshape(mcal_F_rw(:,index_t_origin_2022Q1,i),info_m,info_nvar);
    structpara = [A(:);F(:)];
    LIRF       = IRF_horizons(structpara, info_nvar, info_nlag, info_m, 0:Hhistdecomp_2022Q1);
    B          = F/A;
  

    for ell=1:(Hhistdecomp_2022Q1+1)
        % y(t)' A(t) = x(t)' F(t) + eps(t)'; y(t)' = x(t)' F(t) (A(t))^-1 + eps(t)' A(t)^-1; y(t)' = x(t)' B(t) + u(t)'; u(t)' = y(index_t_origin+t,:)-x(index_t_origin+t,:)*B; u(t)' = eps(t)' A(t)^-1; u(t)'*A = eps(t)'    
        shocks_RC_fixed_para_2022Q1(:,ell,i) = ((y(index_t_origin_2022Q1+ell,:)-x(index_t_origin_2022Q1+ell,:)*B)*A)';
    end
    

    for i_var=1:info_nvar
        for j_shock=1:info_nvar
            for h=0:Hhistdecomp_2022Q1
                switch h
                    case 0
                        xtprime = [1,y(index_t_origin_2022Q1,:),y(index_t_origin_2022Q1-1,:)];
                    case 1
                        xtprime = [1,squeeze(fcast_2022Q1(:,h,i))',y(index_t_origin_2022Q1,:)];
                    otherwise
                        xtprime = [1,squeeze(fcast_2022Q1(:,h,i))',squeeze(fcast_2022Q1(:,h-1,i))'];
                end
                fcast_2022Q1(:,h+1,i) = B'*xtprime';
                tmp=0;
                for ell=0:h
                    tmp= tmp + e(:,i_var)'*LIRF(1+ell*info_nvar:(ell+1)*info_nvar,:)*e(:,j_shock)*e(:,j_shock)'*shocks_RC_fixed_para_2022Q1(:,1+h-ell,i);
                end
                shock_cont_2022Q1(i_var,j_shock,h+1,i) = tmp;
            end
        end
    end

    % historical decomposition 2021Q1
    A          = reshape(mcal_A_rw(:,index_t_origin_2021Q1,i),info_nvar,info_nvar);
    F          = reshape(mcal_F_rw(:,index_t_origin_2021Q1,i),info_m,info_nvar);
    structpara = [A(:);F(:)];
    LIRF       = IRF_horizons(structpara, info_nvar, info_nlag, info_m, 0:Hhistdecomp_2021Q1);
    B          = F/A;
  

    for ell=1:(Hhistdecomp_2021Q1+1)   
        shocks_RC_fixed_para_2021Q1(:,ell,i) = ((y(index_t_origin_2021Q1+ell,:)-x(index_t_origin_2021Q1+ell,:)*B)*A)';
    end
    

    for i_var=1:info_nvar
        for j_shock=1:info_nvar
            for h=0:Hhistdecomp_2021Q1
                switch h
                    case 0
                        xtprime = [1,y(index_t_origin_2021Q1,:),y(index_t_origin_2021Q1-1,:)];
                    case 1
                        xtprime = [1,squeeze(fcast_2021Q1(:,h,i))',y(index_t_origin_2021Q1,:)];
                    otherwise
                        xtprime = [1,squeeze(fcast_2021Q1(:,h,i))',squeeze(fcast_2021Q1(:,h-1,i))'];
                end
                fcast_2021Q1(:,h+1,i) = B'*xtprime';
                tmp=0;
                for ell=0:h
                    tmp= tmp + e(:,i_var)'*LIRF(1+ell*info_nvar:(ell+1)*info_nvar,:)*e(:,j_shock)*e(:,j_shock)'*shocks_RC_fixed_para_2021Q1(:,1+h-ell,i);
                end
                shock_cont_2021Q1(i_var,j_shock,h+1,i) = tmp;
            end
        end
    end

end


% counterfactuals
shocks_RC_cf            = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
count_RC_half_p         = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
count_RC_double_p       = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
count_RC_half_p_Lucas   = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
count_RC_double_p_Lucas = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);

for i=1:1:ndraws
    for t=1:(Hhistdecomp_2022Q1+1)

        A = reshape(mcal_A_rw(:,index_t_origin_2022Q1+t,i),info_nvar,info_nvar);
        F = reshape(mcal_F_rw(:,index_t_origin_2022Q1+t,i),info_m,info_nvar);
        B = F/A;

        % first counterfactual
        A_count_half_p        = A;
        A_count_half_p(2,1)   = A(2,1)/2;
        A_count_double_p      = A;
        A_count_double_p(2,1) = A(2,1)*2;
        % second counterfactual
        A_count_half_p_Lucas   = A;
        A_count_double_p_Lucas = A;

        shocks_RC_cf(:,t,i)    = ((y(index_t_origin_2022Q1+t,:)-x(index_t_origin_2022Q1+t,:)*B)*A)';
        eeps                   = (squeeze(shocks_RC_cf(:,t,i))');


        switch t
            case 1

                eeps_half_p_Lucas = eeps;
                eeps_double_p_Lucas = eeps;
                scale_ffr = 0.75/4;
                eeps_half_p_Lucas(1,1)         = eeps(1,1)  - scale_ffr/median(squeeze(Ltilde_RC(1,3,1,:,ttime==2022.25)));
                eeps_double_p_Lucas(1,1)       = eeps(1,1) + scale_ffr/median(squeeze(Ltilde_RC(1,3,1,:,ttime==2022.25)));
                count_RC_half_p(:,t,i)         = x(index_t_origin_2022Q1+t,:)*F/A_count_half_p + eeps/A_count_half_p;
                count_RC_half_p_Lucas(:,t,i)   = x(index_t_origin_2022Q1+t,:)*F/A_count_half_p_Lucas + eeps_half_p_Lucas/A_count_half_p_Lucas;
                count_RC_double_p(:,t,i)       = x(index_t_origin_2022Q1+t,:)*F/A_count_double_p + eeps/A_count_double_p;
                count_RC_double_p_Lucas(:,t,i) = x(index_t_origin_2022Q1+t,:)*F/A_count_double_p_Lucas + eeps_double_p_Lucas/A_count_double_p_Lucas;

            case 2

                xtcount_half_p                 = [1,count_RC_half_p(:,t-1,i)',y(index_t_origin_2022Q1,:)];
                count_RC_half_p(:,t,i)         = xtcount_half_p*F/A_count_half_p + eeps/A_count_half_p;
                xtcount_double_p               = [1,count_RC_double_p(:,t-1,i)',y(index_t_origin_2022Q1,:)];
                count_RC_double_p(:,t,i)       = xtcount_double_p*F/A_count_double_p + eeps/A_count_double_p;
                xtcount_half_p_Lucas           = [1,count_RC_half_p_Lucas(:,t-1,i)',y(index_t_origin_2022Q1,:)];
                count_RC_half_p_Lucas(:,t,i)   = xtcount_half_p_Lucas*F/A_count_half_p_Lucas + eeps/A_count_half_p_Lucas;
                xtcount_double_p_Lucas         = [1,count_RC_double_p_Lucas(:,t-1,i)',y(index_t_origin_2022Q1,:)];
                count_RC_double_p_Lucas(:,t,i) = xtcount_double_p_Lucas*F/A_count_double_p_Lucas + eeps/A_count_double_p_Lucas;

            otherwise

                xtcount_half_p                 = [1,count_RC_half_p(:,t-1,i)',count_RC_half_p(:,t-2,i)'];
                count_RC_half_p(:,t,i)         = xtcount_half_p*F/A_count_half_p + eeps/A_count_half_p;
                xtcount_double_p               = [1,count_RC_double_p(:,t-1,i)',count_RC_double_p(:,t-2,i)'];
                count_RC_double_p(:,t,i)       = xtcount_double_p*F/A_count_double_p + eeps/A_count_double_p;
                xtcount_half_p_Lucas           = [1,count_RC_half_p_Lucas(:,t-1,i)',count_RC_half_p_Lucas(:,t-2,i)'];
                count_RC_half_p_Lucas(:,t,i)   = xtcount_half_p_Lucas*F/A_count_half_p_Lucas + eeps/A_count_half_p_Lucas;
                xtcount_double_p_Lucas         = [1,count_RC_double_p_Lucas(:,t-1,i)',count_RC_double_p_Lucas(:,t-2,i)'];
                count_RC_double_p_Lucas(:,t,i) = xtcount_double_p_Lucas*F/A_count_double_p_Lucas + eeps/A_count_double_p_Lucas;

        end

    end


end


cq16_2ppsip=4*quantile(count_RC_double_p,0.16,3);
cq50_2ppsip=4*quantile(count_RC_double_p,0.5,3);
cq84_2ppsip=4*quantile(count_RC_double_p,0.84,3);

cq16_2ppsip_ycum=1*quantile(cumsum(count_RC_double_p,2),0.16,3);
cq50_2ppsip_ycum=1*quantile(cumsum(count_RC_double_p,2),0.5,3);
cq84_2ppsip_ycum=1*quantile(cumsum(count_RC_double_p,2),0.84,3);


cq16_2ppsip_Lucas=4*quantile(count_RC_double_p_Lucas,0.16,3);
cq50_2ppsip_Lucas=4*quantile(count_RC_double_p_Lucas,0.5,3);
cq84_2ppsip_Lucas=4*quantile(count_RC_double_p_Lucas,0.84,3);

cq16_2ppsip_Lucas_ycum=1*quantile(cumsum(count_RC_double_p_Lucas,2),0.16,3);
cq50_2ppsip_Lucas_ycum=1*quantile(cumsum(count_RC_double_p_Lucas,2),0.5,3);
cq84_2ppsip_Lucas_ycum=1*quantile(cumsum(count_RC_double_p_Lucas,2),0.84,3);

cq16_0dot5ppsip=4*quantile(count_RC_half_p,0.16,3);
cq50_0dot5ppsip=4*quantile(count_RC_half_p,0.5,3);
cq84_0dot5ppsip=4*quantile(count_RC_half_p,0.84,3);

cq16_0dot5ppsip_ycum=1*quantile(cumsum(count_RC_half_p,2),0.16,3);
cq50_0dot5ppsip_ycum=1*quantile(cumsum(count_RC_half_p,2),0.5,3);
cq84_0dot5ppsip_ycum=1*quantile(cumsum(count_RC_half_p,2),0.84,3);

cq16_0dot5ppsip_Lucas=4*quantile(count_RC_half_p_Lucas,0.16,3);
cq50_0dot5ppsip_Lucas=4*quantile(count_RC_half_p_Lucas,0.5,3);
cq84_0dot5ppsip_Lucas=4*quantile(count_RC_half_p_Lucas,0.84,3);


cq16_0dot5ppsip_Lucas_ycum=1*quantile(cumsum(count_RC_half_p_Lucas,2),0.16,3);
cq50_0dot5ppsip_Lucas_ycum=1*quantile(cumsum(count_RC_half_p_Lucas,2),0.5,3);
cq84_0dot5ppsip_Lucas_ycum=1*quantile(cumsum(count_RC_half_p_Lucas,2),0.84,3);





%% Figure 10



%% Computing objects of interest for R1 + R2 in all periods identification scheme
% load results to get dimensions
cd(current_path);
R2all=load(['results', filesep, ['temp_results_jointpr_','rest118r997ndraws500000_restALL', '.mat']]);

nburn00 = 0.25*size(R2all.mat_L0,4);
R2all.mat_L0 = R2all.mat_L0(:,:,:,nburn00+1:end);
R2all.mat_B  = R2all.mat_B(:,:,nburn00+1:end);

R2all.ndraws  = size(R2all.mat_L0,4);
R2all.A_old   = zeros(R2all.info_nvar);
R2all.Ltilde_RC    = nan(info.hbar+1,R2all.info_nvar,R2all.info_nvar,R2all.ndraws,R2all.info_T);
R2all.CLtilde_RC   = nan(info.hbar+1,R2all.info_nvar,R2all.info_nvar,R2all.ndraws,R2all.info_T);
R2all.shocks_RC    = nan(R2all.info_nvar,R2all.ndraws,R2all.info_T);
R2all.mcal_A_rw    = nan(R2all.info_nvar*R2all.info_nvar,R2all.info_T,R2all.ndraws);
R2all.mcal_F_rw    = nan(R2all.info_k,R2all.info_T,R2all.ndraws);

for i = 1:R2all.ndraws

    % impulse responses and structural parameters
    for t=1:info_T
        R2all.A_old(:,:,t)  = inv(squeeze(R2all.mat_L0(:,:,t,i)))';
        R2all.A = squeeze(R2all.A_old(:,:,t));
        R2all.B = reshape(R2all.mat_B(t,:,i),R2all.info_m,R2all.info_nvar);
        R2all.F=R2all.B*R2all.A;
        R2all.structpara = [R2all.A(:);R2all.F(:)];
        R2all.shocks_RC(:,i,t) = ((R2all.y(t,:)-R2all.x(t,:)*R2all.B)*R2all.A)';
        R2all.LIRF       = IRF_horizons(R2all.structpara, R2all.info_nvar, R2all.info_nlag, R2all.info_m, 0:info.hbar);
        for h=0:info.hbar
            R2all.Ltilde_RC(h+1,:,:,i,t) =  R2all.LIRF(1+h*info_nvar:(h+1)*info_nvar,:);
            for i_shock=1:R2all.info_nvar
                R2all.CLtilde_RC(h+1,[1 2 4],i_shock,i,t)   = sum(R2all.Ltilde_RC(1:h+1,[1 2 4],i_shock,i,t),1);
                R2all.CLtilde_RC(h+1,[3 5],i_shock,i,t)     = R2all.Ltilde_RC(h+1,[3 5],i_shock,i,t);
            end
        end
        R2all.mcal_A_rw(:,t,i) = R2all.A(:);
        R2all.mcal_F_rw(:,t,i) = R2all.F(:);
    end

end




%% Computing objects of interest for constant parameters SVAR
% load results to get dimensions
%/Users/c1jea02/Dropbox/time-varying-SVAR/Replication_package/constant_parameters/figure_1/results
cd ..
resultsconstant_path = pwd;
cd(current_path)

CONSP                      = load([resultsconstant_path,'/constant_parameters/code/results/','results.mat']);
CONSP.t_origin             = 2022;
CONSP.ttime                = 1960.25:0.25:2023.25;
CONSP.index_t_origin       = find(CONSP.ttime==CONSP.t_origin);
CONSP.H                    = 4;
CONSP.info_nlag            = CONSP.nlag;
CONSP.nd                   = size(CONSP.Ltilde,4);
CONSP.info_nvar            = size(CONSP.Ltilde,2);
CONSP.info_m               = size(CONSP.Aplustilde,1);
CONSP.fcast                = nan(CONSP.info_nvar,CONSP.H+1,CONSP.nd);
CONSP.shock_cont           = nan(CONSP.info_nvar,CONSP.info_nvar,CONSP.H+1,CONSP.nd);
CONSP.shocks_RC_fixed_para = nan(CONSP.info_nvar,CONSP.H+1,CONSP.nd);
CONSP.shocks_RC            = nan(CONSP.info_nvar,CONSP.nd,length(CONSP.ttime));


CONSP.y = CONSP.Y;
CONSP.x = CONSP.X;

for i=1:1:CONSP.nd
    CONSP.A = CONSP.A0tilde(:,:,i); 
    CONSP.F = CONSP.Aplustilde(:,:,i); 
    CONSP.structpara = [CONSP.A(:);CONSP.F(:)];
    CONSP.LIRF       = IRF_horizons_CONSP(CONSP.structpara, CONSP.info_nvar, CONSP.info_nlag, CONSP.info_m, 0:5);
    CONSP.B = CONSP.F/CONSP.A;
    for t=1:length(CONSP.ttime)
        CONSP.shocks_RC(:,i,t) = ((CONSP.y(t,:)-CONSP.x(t,:)*CONSP.B)*CONSP.A)';
    end
    for t=1:(CONSP.H+1)
        CONSP.shocks_RC_fixed_para(:,t,i) = ((CONSP.y(CONSP.index_t_origin+t,:)-CONSP.x(CONSP.index_t_origin+t,:)*CONSP.B)*CONSP.A)';
    end
   
    for i_var=1:CONSP.info_nvar
        for j_shock=1:CONSP.info_nvar
            for h=0:CONSP.H
                switch h
                    case 0
                        CONSP.xtprime = [CONSP.y(CONSP.index_t_origin,:),CONSP.y(CONSP.index_t_origin-1,:),CONSP.y(CONSP.index_t_origin-2,:),CONSP.y(CONSP.index_t_origin-3,:),1];
                    case 1
                        CONSP.xtprime = [squeeze(CONSP.fcast(:,h,i))',CONSP.y(CONSP.index_t_origin,:),CONSP.y(CONSP.index_t_origin-1,:),CONSP.y(CONSP.index_t_origin-2,:),1];
                   case 2
                        CONSP.xtprime = [squeeze(CONSP.fcast(:,h,i))',squeeze(CONSP.fcast(:,h-1,i))',CONSP.y(CONSP.index_t_origin,:),CONSP.y(CONSP.index_t_origin-1,:),1];
                   case 3
                        CONSP.xtprime = [squeeze(CONSP.fcast(:,h,i))',squeeze(CONSP.fcast(:,h-1,i))',squeeze(CONSP.fcast(:,h-2,i))',CONSP.y(CONSP.index_t_origin,:),1];
                    otherwise
                        CONSP.xtprime = [squeeze(CONSP.fcast(:,h,i))',squeeze(CONSP.fcast(:,h-1,i))',squeeze(CONSP.fcast(:,h-2,i))',squeeze(CONSP.fcast(:,h-3,i))',1];
                end
                CONSP.fcast(:,h+1,i) = CONSP.B'*CONSP.xtprime';
                tmp=0;
                for ell=0:h
                    tmp= tmp + e(:,i_var)'*CONSP.LIRF(1+ell*CONSP.info_nvar:(ell+1)*CONSP.info_nvar,:)*e(:,j_shock)*e(:,j_shock)'*CONSP.shocks_RC_fixed_para(:,1+h-ell,i);
                end
                CONSP.shock_cont(i_var,j_shock,h+1,i) = tmp;
            end
        end
    end
end


CONSP.Ltildeq50=zeros(size(CONSP.Ltilde,1),size(CONSP.Ltilde,2),size(CONSP.Ltilde,3)); % store IRF quantile 50th
CONSP.Ltildeq16=zeros(size(CONSP.Ltilde,1),size(CONSP.Ltilde,2),size(CONSP.Ltilde,3)); % store IRF quantile 16th
CONSP.Ltildeq84=zeros(size(CONSP.Ltilde,1),size(CONSP.Ltilde,2),size(CONSP.Ltilde,3)); % store IRF quantile 84th

CONSP.CLtildeq50=zeros(size(CONSP.CLtilde,1),size(CONSP.CLtilde,2),size(CONSP.CLtilde,3)); % store IRF quantile 50th
CONSP.CLtildeq16=zeros(size(CONSP.CLtilde,1),size(CONSP.CLtilde,2),size(CONSP.CLtilde,3)); % store IRF quantile 16th
CONSP.CLtildeq84=zeros(size(CONSP.CLtilde,1),size(CONSP.CLtilde,2),size(CONSP.CLtilde,3)); % store IRF quantile 84th

CONSP.numdate_to_plot=2022.5;
CONSP.scale_shock = 4*squeeze(median(squeeze(CLtilde_RC(:,3,1,:,ttime==CONSP.numdate_to_plot)),2));
for ii=1:size(CONSP.Ltilde,1)
    for jj=1:size(CONSP.Ltilde,2)
        for kk=1:size(CONSP.Ltilde,3)
            
            CONSP.Ltildeq50(ii,jj,kk) = quantile(CONSP.Ltilde(ii,jj,kk,:)*CONSP.scale_shock(1)/median(CONSP.Ltilde(1,3,1,:)),0.5);
            CONSP.Ltildeq16(ii,jj,kk) = quantile(CONSP.Ltilde(ii,jj,kk,:)*CONSP.scale_shock(1)/median(CONSP.Ltilde(1,3,1,:)),0.16);
            CONSP.Ltildeq84(ii,jj,kk) = quantile(CONSP.Ltilde(ii,jj,kk,:)*CONSP.scale_shock(1)/median(CONSP.Ltilde(1,3,1,:)),0.84);
            
            
            CONSP.CLtildeq50(ii,jj,kk) = quantile(CONSP.CLtilde(ii,jj,kk,:)*CONSP.scale_shock(1)/median(CONSP.CLtilde(1,3,1,:)),0.5);
            CONSP.CLtildeq16(ii,jj,kk) = quantile(CONSP.CLtilde(ii,jj,kk,:)*CONSP.scale_shock(1)/median(CONSP.CLtilde(1,3,1,:)),0.16);
            CONSP.CLtildeq84(ii,jj,kk) = quantile(CONSP.CLtilde(ii,jj,kk,:)*CONSP.scale_shock(1)/median(CONSP.CLtilde(1,3,1,:)),0.84);
            
            
        end
    end
end



%% Figure III.1
gcafontsize=14;
hFig=figure('name','shocks');
set(hFig, 'Position', [0 20 650 350])
a16_conj = quantile(squeeze(shocks_RC(1,:,210:215))',0.16,2)';
b84_conj = quantile(squeeze(shocks_RC(1,:,210:215))',0.84,2)';
[~,~]=jbfill(ttime(210:215),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.75);
hold on
plot(ttime(210:215),median(squeeze(shocks_RC(1,:,210:215)),1),'color',rgb('darkgreen'),'linewidth',1,'LineStyle','-')
hold on
hline(0,':k')
hold on
hold on
xline(ttime(ttime==2022.5),'--r','2022Q3 (Contractionary)','fontsize',14)
xlim([ttime(210) ttime(end)])
xticks([ttime(210) ttime(211) ttime(212) ttime(213) ttime(214) ttime(215)])
grid on
xticklabels({'2022Q1','2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
hold on
% recessionplot
 hold on
%title('Monetary Policy Shocks')
box off
grid on
set(gca,'Fontsize',gcafontsize)


print([savepath, filesep, 'figure_III_1.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_III_1.png'],'-dpng');




% Figure III.6 panel a
gcafontsize=14;
hFig=figure('name','Figure III.6 panel a');
set(hFig, 'Position', [0 20 450 300])
ttime_start=206;
a16_conj = quantile(squeeze(shocks_RC(1,1:end,ttime_start:end))',0.16,2)';
b84_conj = quantile(squeeze(shocks_RC(1,1:end,ttime_start:end))',0.84,2)';
R2all.a16_conj = quantile(squeeze(R2all.shocks_RC(1,1:end,ttime_start:end))',0.16,2)';
R2all.b84_conj = quantile(squeeze(R2all.shocks_RC(1,1:end,ttime_start:end))',0.84,2)';
[~,~]=jbfill(ttime(ttime_start:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.75);
hold on
[~,~]=jbfill(ttime(ttime_start:end),R2all.a16_conj,R2all.b84_conj,rgb('lightgray'),rgb('lightgray'),0,0.75);
hold on
plot(ttime(ttime_start:end),median(squeeze(shocks_RC(1,:,ttime_start:end)),1),'color',rgb('darkgreen'),'linewidth',1,'LineStyle','-')
hold on
plot(ttime(ttime_start:end),median(squeeze(R2all.shocks_RC(1,:,ttime_start:end)),1),'color',rgb('darkgray'),'linewidth',1,'LineStyle','-')
hold on
hline(0,':k')
hold on
hold on
xlim([ttime(ttime_start) ttime(end)])
hold on
recessionplot
hold on
box off
grid on
xticks([ttime(ttime_start) ttime(ttime_start+2) ttime(ttime_start+4)  ttime(ttime_start+6) ttime(ttime_start+8)])
grid on
xticklabels({'2021Q1','2021Q3','2022Q1','2022Q3','2023Q1'})
set(gca,'Fontsize',gcafontsize)

print([savepath, filesep, 'figure_III_6a.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_III_6a.png'],'-dpng');




% Figure III.6 panel b
hFig=figure('name','Figure III.6 panel b');
set(hFig, 'Position', [0 20 450 300])
histogram(squeeze(shocks_RC(1,:,ttime_start+3)),'facecolor',rgb('lightgreen'),'normalization','probability')
hold on
histogram(squeeze(R2all.shocks_RC(1,:,ttime_start+3)),'facecolor',rgb('lightgray'),'normalization','probability')
box off
xline(0,':k','linewidth',2)
set(gca,'Fontsize',gcafontsize)
grid on

print([savepath, filesep, 'figure_III_6b.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_III_6b.png'],'-dpng');

% Figure III.6 panel c
hFig=figure('name','Figure III.6 panel c');
set(hFig, 'Position', [0 20 450 300])
numdate_to_plot=2021.25;
h=histogram(-1.*squeeze(CLtilde_RC(1,1,1,:,ttime==numdate_to_plot)),'facecolor',rgb('lightgreen'),'normalization','probability');
hold on
histogram(-1*squeeze(R2all.CLtilde_RC(1,1,1,:,ttime==numdate_to_plot)),'facecolor',rgb('lightgray'),'normalization','probability','BinEdges',h.BinEdges)
box off
grid on
xline(0,':k','linewidth',2)
ylim([0 0.1])
xlim([-3,3])
set(gca,'Fontsize',14)

print([savepath, filesep, 'figure_III_6c.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_III_6c.png'],'-dpng');




% figure III.7
dates2plot=2022.25:0.25:2023.25;
hFig = figure('name','Figure III_7');
set(hFig, 'Position', [0 250 1000 250])
subplot(1,3,1)
tmp=[4*mean(squeeze(CONSP.shock_cont(3,1,:,:)),2),4*mean(squeeze(sum(squeeze(CONSP.shock_cont(3,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
hold on
leg=legend([h(1) h(2)],'Monetary Policy Shock','Non-Monetary Policy Shocks','NumColumns', 1);
set(leg,'Location','northwest',...
    'Orientation','vertical',...
    'NumColumns',1,...
    'FontSize',10);
set(gca,'Fontsize',12)
xtickangle(90)
yticks([0 1 2 3 4])
ylim([0,4])
ylabel('% (annualized)')
title('Federal Funds Rate')
legend boxoff
set(gcf, 'PaperPositionMode', 'auto');

subplot(1,3,2)
tmp=[4*mean(squeeze(CONSP.shock_cont(1,1,:,:)),2),4*mean(squeeze(sum(squeeze(CONSP.shock_cont(1,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-3 -2 -1 0 1 2 3 4])
ylim([-3,4])
ylabel('% (log-difference annualized)')
title('Output Growth')

subplot(1,3,3)
hold on
hold on
tmp=[4*mean(squeeze(CONSP.shock_cont(2,1,:,:)),2),4*mean(squeeze(sum(squeeze(CONSP.shock_cont(2,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
hold on
hold on
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
title('Core Inflation')
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-1.5 -0.5 0 0.5 1 1.5 2])
ylim([-1.5,2])
ylabel('% (log-difference annualized)')



print([savepath, filesep, 'Figure_III_7a.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'Figure_III_7a.png'],'-dpng');



%Figure III_7_b
dates2plot=2022.25:0.25:2023.25;
hFig = figure('name','Figure III_7b');
set(hFig, 'Position', [0 250 1000 250])

subplot(1,3,1)
tmp=[4*mean(squeeze(shock_cont_2022Q1(3,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2022Q1(3,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
hold on
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
hold on
hold on
leg=legend([h(1) h(2)],'Monetary Policy Shock','Non-Monetary Policy Shocks', 'NumColumns', 1);
set(leg,'Location','north',...
    'Orientation','vertical',...
    'NumColumns',1,...
    'FontSize',10);
set(gca,'Fontsize',12)
xtickangle(90)
yticks([0 1 2 3 4])
ylim([0,4])
ylabel('% (annualized)')
title('Federal Funds Rate')
legend boxoff
set(gcf, 'PaperPositionMode', 'auto');

subplot(1,3,2)
tmp=[4*mean(squeeze(shock_cont_2022Q1(1,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2022Q1(1,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-3 -2 -1 0 1 2 3 4])
ylim([-3,4])
ylabel('% (log-difference annualized)')
title('Output Growth')

subplot(1,3,3)
tmp=[4*mean(squeeze(shock_cont_2022Q1(2,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2022Q1(2,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
title('Core Inflation')
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-1.5 -0.5 0 0.5 1 1.5 2])
ylim([-1.5,2])
ylabel('% (log-difference annualized)')



print([savepath, filesep, 'Figure_III_7b.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'Figure_III_7b.png'],'-dpng');




% Figure III_8
ftsizeaxis=14;
ftsizexlabel=12;
ftsizetitle=14;
ftlinewidth = 1.0;
medianwidth=2.0;

hFig = figure('name','IRF 2022Q3');
set(hFig, 'Position', [0 250 1000 250])
addpath('figures_helpfunctions/')
bands_color='gray';
CONSP.HIRF = size(CONSP.Ltildeq16,1);
subplot(1,3,1)
x = (0:1:(CONSP.HIRF-1));
hold on
[~,~]=jbfill(x,4*quantile(squeeze(CLtilde_RC(:,3,1,:,ttime==CONSP.numdate_to_plot)),0.16,2)',4*quantile(squeeze(CLtilde_RC(:,3,1,:,ttime==CONSP.numdate_to_plot)),0.84,2)',rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
plot(x,4*quantile(squeeze(CLtilde_RC(:,3,1,:,ttime==CONSP.numdate_to_plot)),0.5,2)','color',rgb('darkgreen'),'LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.Ltildeq16(1:CONSP.HIRF,3,1)),'color',rgb('black'),'LineStyle',':','LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.Ltildeq84(1:CONSP.HIRF,3,1)),'color',rgb('black'),'LineStyle',':','LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.Ltildeq50(1:CONSP.HIRF,3,1)),'color',rgb('black'),'LineWidth',medianwidth)
hold on
hline(0,'-r')
hold on
set(gca,'LineWidth',ftlinewidth )
set(gca,'YTick',[-0.5 -0.25 0 0.25 0.5 0.75 1])
axis([0 (CONSP.HIRF-1) -0.5 1])
xlabel('Quarters','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Federal Funds Rate','FontSize',ftsizetitle)
box off
grid on


subplot(1,3,2)
[~,~]=jbfill(x,quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==CONSP.numdate_to_plot)),0.16,2)',quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==CONSP.numdate_to_plot)),0.84,2)',rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
plot(x,quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==CONSP.numdate_to_plot)),0.5,2)','color',rgb('darkgreen'),'LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.CLtildeq16(1:CONSP.HIRF,1,1)),'color',rgb('black'),'LineStyle',':','LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.CLtildeq84(1:CONSP.HIRF,1,1)),'color',rgb('black'),'LineStyle',':','LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.CLtildeq50(1:CONSP.HIRF,1,1)),'color',rgb('black'),'LineWidth',medianwidth)
hold on
hold on
hline(0,'-r')
hold on
set(gca,'LineWidth',ftlinewidth )
set(gca,'YTick',[-4 -3 -2 -1 0 1 2 3 4])
axis([0 (CONSP.HIRF-1) -4 4])
xlabel('Quarters','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Output','FontSize',ftsizetitle)
box off
grid on


subplot(1,3,3)
[~,~]=jbfill(x,quantile(4*squeeze(Ltilde_RC(:,2,1,:,ttime==CONSP.numdate_to_plot)),0.16,2)',quantile(4*squeeze(Ltilde_RC(:,2,1,:,ttime==CONSP.numdate_to_plot)),0.84,2)',rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
hold on
plot(x,quantile(4*squeeze(Ltilde_RC(:,2,1,:,ttime==CONSP.numdate_to_plot)),0.5,2)','color',rgb('darkgreen'),'LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.Ltildeq16(1:CONSP.HIRF,2,1)),'color',rgb('black'),'LineStyle',':','LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.Ltildeq84(1:CONSP.HIRF,2,1)),'color',rgb('black'),'LineStyle',':','LineWidth',medianwidth)
hold on
plot(x,squeeze(CONSP.Ltildeq50(1:CONSP.HIRF,2,1)),'color',rgb('black'),'LineWidth',medianwidth)
hold on
hold on
hline(0,'-r')
hold on
set(gca,'LineWidth',ftlinewidth )
xlabel('Quarters','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Core Inflation','FontSize',ftsizetitle)
box off
grid on

set(gcf, 'PaperPositionMode', 'auto');



print([savepath, filesep, 'Figure_III_8.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'Figure_III_8.png'],'-dpng');




% Figure III.9
dates2plot=2022.25:0.25:2023.25;
hFig = figure('name','Figure 11');
set(hFig, 'Position', [0 250 1000 250])
subplot(1,3,1)
tmp=[4*mean(squeeze(CONSP.shock_cont(3,1,:,:)),2),4*mean(squeeze(sum(squeeze(CONSP.shock_cont(3,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
leg=legend([h(1) h(2)],'Monetary Policy Shock','Non-Monetary Policy Shocks','NumColumns', 1);
set(leg,'Location','northwest',...
    'Orientation','vertical',...
    'NumColumns',1,...
    'FontSize',10);
set(gca,'Fontsize',12)
xtickangle(90)
yticks([0 1 2 3 4])
ylim([0,4])
ylabel('% (annualized)')
title('Federal Funds Rate')
legend boxoff
set(gcf, 'PaperPositionMode', 'auto');

subplot(1,3,2)
tmp=[4*mean(squeeze(CONSP.shock_cont(4,1,:,:)),2),4*mean(squeeze(sum(squeeze(CONSP.shock_cont(4,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-12 -9 -6 -3 0])
ylim([-12,0])
ylabel('% (log-difference annualized)')
title('Money Growth')

subplot(1,3,3)

tmp=[4*mean(squeeze(CONSP.shock_cont(5,1,:,:)),2),4*mean(squeeze(sum(squeeze(CONSP.shock_cont(5,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-0.4 -0.2 0 0.2 0.4])
ylim([-0.4,0.4])
ylabel('% (annualized)')
title(' Credit Spreads')

print([savepath, filesep, 'figure_III_9.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_III_9.png'],'-dpng');







% Figure IV.1
gcafontsize=12;
hFig=figure('name','Robust');
set(hFig, 'Position', [0 20 950 300])

subplot(1,3,1)
a16_conj=cq16_2ppsip_Lucas(3,:);
b84_conj=cq84_2ppsip_Lucas(3,:);
[p1,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
a16_conj=cq16_0dot5ppsip_Lucas(3,:);
b84_conj=cq84_0dot5ppsip_Lucas(3,:);
[p2,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('khaki'),rgb('khaki'),0,0.55);
hold on
p3=plot(ttime(index_t_origin_2022Q1+1:end),y(end-4:end,3)*4,'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip_Lucas(3,:),'linewidth',2,'color',rgb('green'))
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip_Lucas(3,:),'linewidth',2,'color',rgb('yellow'))
axis('tight')
leg=legend([p1 p2 p3],'Hawkish Fed','Dovish Fed','Data');
set(leg,'location','best','fontsize',10);
legend boxoff
hold on
hline(0,':r')
box off
ylabel('%')
grid on
title('Federal Funds Rate')
xticks([2022.25 2022.5 2022.75 2023 2023.25])
yticks([0 2 4 6 8 10])
ylim([-0.5,10])
xlim([2022.25,2023.25])
xticklabels({'22Q2','22Q3','22Q4','23Q1','23Q2'})
set(gca,'Fontsize',gcafontsize)
xtickangle(0)


subplot(1,3,2)
a16_conj=cq16_2ppsip_Lucas_ycum(1,:);
b84_conj=cq84_2ppsip_Lucas_ycum(1,:);
jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
a16_conj=cq16_0dot5ppsip_Lucas_ycum(1,:);
b84_conj=cq84_0dot5ppsip_Lucas_ycum(1,:);
jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('khaki'),rgb('khaki'),0,0.5);
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cumsum(y(end-4:end,1)*1),'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'))
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip_Lucas_ycum(1,:),'linewidth',2,'color',rgb('green'))
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip_Lucas_ycum(1,:),'linewidth',2,'color',rgb('yellow'))
hold on
hline(0,':r')
box off
title('Output')
ylabel('% (log-change)')
grid on
xticks([2022.25 2022.5 2022.75 2023 2023.25])
yticks([-4 -2 0 2 4 6])
ylim([-4,6])
xlim([2022.25,2023.25])
xticklabels({'22Q2','22Q3','22Q4','23Q1','23Q2'})
set(gca,'Fontsize',gcafontsize)
xtickangle(0)

subplot(1,3,3)
a16_conj=cq16_2ppsip_Lucas(2,:);
b84_conj=cq84_2ppsip_Lucas(2,:);
jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
a16_conj=cq16_0dot5ppsip_Lucas(2,:);
b84_conj=cq84_0dot5ppsip_Lucas(2,:);
jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('khaki'),rgb('khaki'),0,0.55);
hold on
plot(ttime(index_t_origin_2022Q1+1:end),y(end-4:end,2)*4,'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip_Lucas(2,:),'linewidth',2,'color',rgb('green'))
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip_Lucas(2,:),'linewidth',2,'color',rgb('yellow'))
axis('tight')
hold on
hline(0,':r')
box off
ylabel('%(log-difference annualized)')
grid on
title('Core Inflation')
xticks([2022.25 2022.5 2022.75 2023 2023.25])
yticks([0 2 4 6 8])
ylim([0,8])
xlim([2022.25,2023.25])
xticklabels({'22Q2','22Q3','22Q4','23Q1','23Q2'})
set(gca,'Fontsize',gcafontsize)
xtickangle(0)


 print([savepath, filesep, 'figure_IV_1.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_IV_1.png'],'-dpng');
