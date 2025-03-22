%% Housekeeping
clc; clear variables; close all;


%% Load results 
current_path = pwd;
addpath(genpath('figures_helpfunctions'));
cd ..
addpath(genpath('dcd_helpfunctions'));
results_path = pwd;
cd(current_path);

specname = 'rest115r3024ndraws200000';%'';
savepath = [pwd, filesep, 'results', filesep, specname];

load(['results', filesep, ['temp_results_jointpr_',specname, '.mat']]);

mkdir(['results', filesep, specname]);
savepath = [pwd, filesep, 'results', filesep, specname];
%% burn 1/4 of mcmc chain
max_draws = size(mat_Q,4);
nburn00 = max_draws*0.25;

mat_B = mat_B(:,:,nburn00+1:max_draws);
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
ttime                        = 1969.75:0.25:2023.25;
index_t_origin_2021Q1        = find(ttime==t_origin_2021Q1);
Hhistdecomp_2021Q1           = 8;
fcast_2021Q1                 = nan(info_nvar,Hhistdecomp_2021Q1+1,ndraws);
shock_cont_2021Q1            = nan(info_nvar,info_nvar,Hhistdecomp_2021Q1+1,ndraws);
shocks_RC_fixed_para_2021Q1  = nan(info_nvar,Hhistdecomp_2021Q1+1,ndraws);



for i = 1:ndraws

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
                CLtilde_RC(h+1,[1 2],i_shock,i,t)   = sum(Ltilde_RC(1:h+1,[1 2],i_shock,i,t),1);
                CLtilde_RC(h+1,[3],i_shock,i,t)     = Ltilde_RC(h+1,[3],i_shock,i,t);
            end
        end
        mcal_A_rw(:,t,i) = A(:);
        mcal_F_rw(:,t,i) = F(:);
        
        
        
    end


 

end



%% Appendix


%% Figure V1 
% ru



%% Figure V2
ftsizeaxis=14;
ftsizexlabel=12;
ftsizetitle=14;
ftlinewidth = 1.0;
medianwidth=2.0;

hFig = figure('name','IRF 2022Q1');
set(hFig, 'Position', [0 250 750 300])

numdate_to_plot_2022Q1 = 2022;
numdate_to_plot_1975Q1 = 1975;

bands_color='gray';
HIRF = size(CLtilde_RC,1);

subplot(1,2,1)
[h1,~]=jbfill(x,quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==numdate_to_plot_2022Q1)),0.16,2)',quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==numdate_to_plot_2022Q1)),0.84,2)',rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
plot(x,quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==numdate_to_plot_2022Q1)),0.5,2)','color',rgb('darkgreen'),'LineWidth',medianwidth)
hold on
hold on
[h2,~]=jbfill(x,quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==numdate_to_plot_1975Q1)),0.16,2)',quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==numdate_to_plot_1975Q1)),0.84,2)',rgb('khaki'),rgb('khaki'),0,0.5);
hold on
plot(x,quantile(squeeze(CLtilde_RC(:,1,1,:,ttime==numdate_to_plot_1975Q1)),0.5,2)','color',rgb('yellow'),'LineWidth',medianwidth)
hold on
hline(0,'-r')
hold on
set(gca,'LineWidth',ftlinewidth )
% set(gca,'YTick',[-3 -2 -1 0 1 2 3])
% axis([0 (HIRF-1) -3 3])
xlabel('Quarters','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Output','FontSize',ftsizetitle)
box off
grid on
leg=legend([h1 h2],'2022Q1','1975Q1')
set(leg,'location','northwest')
legend box off

subplot(1,2,2)
[~,~]=jbfill(x,quantile(squeeze(CLtilde_RC(:,2,1,:,ttime==numdate_to_plot_2022Q1)),0.16,2)',quantile(squeeze(CLtilde_RC(:,2,1,:,ttime==numdate_to_plot_2022Q1)),0.84,2)',rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
plot(x,quantile(squeeze(CLtilde_RC(:,2,1,:,ttime==numdate_to_plot_2022Q1)),0.5,2)','color',rgb('darkgreen'),'LineWidth',medianwidth)
hold on
hold on
[~,~]=jbfill(x,quantile(squeeze(CLtilde_RC(:,2,1,:,ttime==numdate_to_plot_1975Q1)),0.16,2)',quantile(squeeze(CLtilde_RC(:,2,1,:,ttime==numdate_to_plot_1975Q1)),0.84,2)',rgb('khaki'),rgb('khaki'),0,0.5);
hold on
plot(x,quantile(squeeze(CLtilde_RC(:,2,1,:,ttime==numdate_to_plot_1975Q1)),0.5,2)','color',rgb('yellow'),'LineWidth',medianwidth)
hold on
hline(0,'-r')
hold on
set(gca,'LineWidth',ftlinewidth )
xlabel('Quarters','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Price Level','FontSize',ftsizetitle)
box off
grid on

set(gcf, 'PaperPositionMode', 'auto');



print([savepath, filesep, 'figure_V2_.eps'],'-depsc');
print([savepath, filesep, 'figure_V2.png'],'-dpng');
