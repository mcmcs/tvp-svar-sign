function [output]=pull_figure_5()
%% Housekeeping
clc; clear variables; close all;


%% Load results 
current_path = pwd;
addpath(genpath('figures_helpfunctions'));
addpath(genpath('tvsvar_helpfunctions'));
addpath(genpath('var_helpfunctions'));
results_path = pwd;
cd(current_path);

% load(['results', filesep, 'temp_results_rest4_zlbcrt20K_com1.mat']);
% load(['results', filesep, 'temp_results_jointpr700.mat']);
% load(['results', filesep, 'temp_results_jointpr_rest2r700.mat']);

specname = 'rest115r3024ndraws1000000';%'';

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
mat_logL = mat_logL(nburn00+1:max_draws,1);


% % find highest likelihood
% [max_logL,ind_max_logL]=max(mat_logL);
% 
% mean_value = mean(mat_logL);          % Compute the median of the vector
% [~, ind_mean_logL] = min(abs(mat_logL - mean_value)); % Find the index closest to the median
% 
% ind_logL= ind_mean_logL;

L0_old = nan(size(mat_L0,1),size(mat_L0,2),size(mat_L0,3));




%% Computing objects of interest for baseline identification scheme
%  1.1. Get the structural parameters
%  1.2. Compute the systematic component of monetary policy
%  1.3. Compute the historical decomposition
%ndraws  = size(mat_L0,4);
n_Q     = 50000;

A_old   = zeros(info_nvar);
e       = eye(info_nvar);

r2.Ssigma    = zeros(info_nvar,info_nvar,n_Q);
info.hbar    = 4; % max IRF horizon
Ltilde_RC    = nan(info.hbar+1,info_nvar,info_nvar,n_Q,info_T);
CLtilde_RC   = nan(info.hbar+1,info_nvar,info_nvar,n_Q,info_T);
shocks_RC    = nan(info_nvar,n_Q,info_T);
mcal_A_rw    = nan(info_nvar*info_nvar,info_T,n_Q);
mcal_F_rw    = nan(info_k,info_T,n_Q);


ttime                        = 1969.75:0.25:2023.25;



% Historical decomposition starting in 2022Q1
t_origin_2022Q1               = 2022;
index_t_origin_2022Q1        = find(ttime==t_origin_2022Q1);

Hhistdecomp_2022Q1           = 4;

function_restrictions = @restriction_baseline_taylor_1s_gt0_gt1;
is_restrict_B = true; %true if restriction includes B  

for i = 1:n_Q

    disp([num2str(i), ' / ', num2str(n_Q)]);


    
    % impulse responses and structural parameters
    for t=index_t_origin_2022Q1+1:1:info_T

        % reduced-form para with highest likelihood
        %Bhat reshape(mat_B(t,:,ind_logL),info_m,info_nvar);
        Bhat = median(mat_B(t,:,:),3); % minchul fixed
        ddeltahat = median(mat_ddelta(t,:,:),3);
        ggammahat = median(mat_ggamma(t,:,:),3);
  
%         disp('code gets stuck here')
% t
        Q_old = draw_Q_condi(Bhat, ddeltahat, ggammahat, function_restrictions, info_nvar, t);


         % IRF on impact

         Dt_old    = diag(exp(ddeltahat/2));
         Ct_old    = ggammatoC(ggammahat);
         DCDhalf_old =  chol(Dt_old*Ct_old*Dt_old, 'lower');

         L0_old(:,:,t) = DCDhalf_old * Q_old;
  

        A_old(:,:,t)  = inv(squeeze(L0_old(:,:,t)))';
        A = squeeze(A_old(:,:,t));
        B = reshape(Bhat,info_m,info_nvar);
        F=B*A;
        mcal_A_rw(:,t,i) = A(:);
        mcal_F_rw(:,t,i) = F(:);
    end
end



ndraws=n_Q;
% counterfactuals
shocks_RC_cf            = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
count_RC_half_p         = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);
count_RC_double_p       = nan(info_nvar,Hhistdecomp_2022Q1+1,ndraws);

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

        shocks_RC_cf(:,t,i)    = ((y(index_t_origin_2022Q1+t,:)-x(index_t_origin_2022Q1+t,:)*B)*A)';
        eeps                   = (squeeze(shocks_RC_cf(:,t,i))');


        switch t
            case 1

     
        
                count_RC_half_p(:,t,i)         = x(index_t_origin_2022Q1+t,:)*F/A_count_half_p + eeps/A_count_half_p;
                count_RC_double_p(:,t,i)       = x(index_t_origin_2022Q1+t,:)*F/A_count_double_p + eeps/A_count_double_p;
       
            case 2

                xtcount_half_p                 = [1,count_RC_half_p(:,t-1,i)',y(index_t_origin_2022Q1,:)];
                count_RC_half_p(:,t,i)         = xtcount_half_p*F/A_count_half_p + eeps/A_count_half_p;
                xtcount_double_p               = [1,count_RC_double_p(:,t-1,i)',y(index_t_origin_2022Q1,:)];
                count_RC_double_p(:,t,i)       = xtcount_double_p*F/A_count_double_p + eeps/A_count_double_p;
      
            otherwise

                xtcount_half_p                 = [1,count_RC_half_p(:,t-1,i)',count_RC_half_p(:,t-2,i)'];
                count_RC_half_p(:,t,i)         = xtcount_half_p*F/A_count_half_p + eeps/A_count_half_p;
                xtcount_double_p               = [1,count_RC_double_p(:,t-1,i)',count_RC_double_p(:,t-2,i)'];
                count_RC_double_p(:,t,i)       = xtcount_double_p*F/A_count_double_p + eeps/A_count_double_p;
        
        end

    end


end

sel_count = intersect(find(count_RC_half_p(3,1,:)>0),find(count_RC_half_p(3,2,:)>0));
sel_count = intersect(sel_count,find(count_RC_half_p(3,3,:)>0));
sel_count = intersect(sel_count,find(count_RC_half_p(3,4,:)>0));
sel_count = intersect(sel_count,find(count_RC_half_p(3,5,:)>0));




cq16_2ppsip=4*quantile(count_RC_double_p,0.16,3);
cq50_2ppsip=4*quantile(count_RC_double_p,0.5,3);
cq84_2ppsip=4*quantile(count_RC_double_p,0.84,3);
cmin_2ppsip=4*min(count_RC_double_p,[],3);
cmax_2ppsip=4*max(count_RC_double_p,[],3);

cq16_2ppsip_ycum=1*quantile(cumsum(count_RC_double_p,2),0.16,3);
cq50_2ppsip_ycum=1*quantile(cumsum(count_RC_double_p,2),0.5,3);
cq84_2ppsip_ycum=1*quantile(cumsum(count_RC_double_p,2),0.84,3);
cmin_2ppsip_ycum=1*min(cumsum(count_RC_double_p,2),[],3);
cmax_2ppsip_ycum=1*max(cumsum(count_RC_double_p,2),[],3);


cq16_0dot5ppsip=4*quantile(count_RC_half_p(:,:,sel_count),0.16,3);
cq50_0dot5ppsip=4*quantile(count_RC_half_p(:,:,sel_count),0.5,3);
cq84_0dot5ppsip=4*quantile(count_RC_half_p(:,:,sel_count),0.84,3);
cmin_0dot5ppsip=4*min(count_RC_half_p(:,:,sel_count),[],3);
cmax_0dot5ppsip=4*max(count_RC_half_p(:,:,sel_count),[],3);


cq16_0dot5ppsip_ycum=1*quantile(cumsum(count_RC_half_p(:,:,sel_count),2),0.16,3);
cq50_0dot5ppsip_ycum=1*quantile(cumsum(count_RC_half_p(:,:,sel_count),2),0.5,3);
cq84_0dot5ppsip_ycum=1*quantile(cumsum(count_RC_half_p(:,:,sel_count),2),0.84,3);
cmin_0dot5ppsip_ycum=1*min(cumsum(count_RC_half_p(:,:,sel_count),2),[],3);
cmax_0dot5ppsip_ycum=1*max(cumsum(count_RC_half_p(:,:,sel_count),2),[],3);




%% Figure 5
gcafontsize=12;
hFig=figure('name','counterfactuals');
set(hFig, 'Position', [0 20 950 300])
save('figure_5_identified_set.mat',...
    'cmin_2ppsip','cmax_2ppsip','cmin_2ppsip_ycum','cmax_2ppsip_ycum',...
    'cmin_0dot5ppsip','cmax_0dot5ppsip','cmin_0dot5ppsip_ycum','cmax_0dot5ppsip_ycum',...
  'ttime','index_t_origin_2022Q1','y','cq50_0dot5ppsip_ycum','cq50_0dot5ppsip','cq50_2ppsip_ycum','cq50_2ppsip')
subplot(1,3,1)
amin_conj=cmin_2ppsip(3,:);
bmax_conj=cmax_2ppsip(3,:);
[p1,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),amin_conj,bmax_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
amin_conj=cmin_0dot5ppsip(3,:);
bmax_conj=cmax_0dot5ppsip(3,:);
[p2,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),amin_conj,bmax_conj,rgb('khaki'),rgb('khaki'),0,0.55);
hold on
p3=plot(ttime(index_t_origin_2022Q1+1:end),y(end-4:end,3)*4,'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
%plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip(3,:),'linewidth',2,'color',rgb('green'))
hold on
%plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip(3,:),'linewidth',2,'color',rgb('yellow'))
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
amin_conj=cmin_2ppsip_ycum(1,:);
bmax_conj=cmax_2ppsip_ycum(1,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),amin_conj,bmax_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
amin_conj=cmin_0dot5ppsip_ycum(1,:);
bmax_conj=cmax_0dot5ppsip_ycum(1,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),amin_conj,bmax_conj,rgb('khaki'),rgb('khaki'),0,0.5);
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cumsum(y(end-4:end,1)*1),'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
%plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip_ycum(1,:),'linewidth',2,'color',rgb('green'))
hold on
%plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip_ycum(1,:),'linewidth',2,'color',rgb('yellow'))
hold on
hline(0,':r')
box off
title('Output')
ylabel('% (log-change)')
grid on
xticks([2022.25 2022.5 2022.75 2023 2023.25])
%yticks([-4 -2 0 2 4 6])
ylim([-2,7])
xlim([2022.25,2023.25])
xticklabels({'22Q2','22Q3','22Q4','23Q1','23Q2'})
set(gca,'Fontsize',gcafontsize)
xtickangle(0)
subplot(1,3,3)
amin_conj=cmin_2ppsip(2,:);
bmax_conj=cmax_2ppsip(2,:);

[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),amin_conj,bmax_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
amin_conj=cmin_0dot5ppsip(2,:);
bmax_conj=cmax_0dot5ppsip(2,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),amin_conj,bmax_conj,rgb('khaki'),rgb('khaki'),0,0.55);
hold on
plot(ttime(index_t_origin_2022Q1+1:end),y(end-4:end,2)*4,'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
%plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip(2,:),'linewidth',2,'color',rgb('green'));
hold on
%plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip(2,:),'linewidth',2,'color',rgb('yellow'));
axis('tight')
hold on
hline(0,':r')
box off
ylabel('%(log-difference annualized)')
grid on
title('Core Inflation')
xticks([2022.25 2022.5 2022.75 2023 2023.25])
%yticks([0 2 4 6 8 10])
ylim([0,18])
xlim([2022.25,2023.25])
xticklabels({'22Q2','22Q3','22Q4','23Q1','23Q2'})
set(gca,'Fontsize',gcafontsize)
xtickangle(0)


% print('results/figure_5.eps','-depsc');
% pause(0.1);
% print('results/figure_5.png','-dpng');


print([savepath, filesep, 'figure_5.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_5.png'],'-dpng');

output=1;

end