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
info.hbar    = 4; % max IRF horizon
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
        LIRF       = IRF_horizons(structpara, info_nvar, info_nlag, info_m, 0:Hhistdecomp_2022Q1);
        for h=0:Hhistdecomp_2022Q1
            
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







%% Figure (1) Systematic Component
ttime=1969.75:0.25:2023.25;
end_date=2023.25;
ttime_m = [1979.75:0.25:1982.75,2009:0.25:2015.5,2020.25:0.25:2021.75];
index_m=ismember(ttime,ttime_m);
close all
gcafontsize=12;
hFig=figure('name','Systematic Component');
set(hFig, 'Position', [20 20 700 500])
subplot(2,2,1)
hold on
plot(ttime,ttime*0,':r')
hold on
a16_conj = quantile(-squeeze(mcal_A_rw(1,:,:))./squeeze(mcal_A_rw(3,:,:)),0.16,2)';
a16_conj(index_m)=nan;
b84_conj = quantile(-squeeze(mcal_A_rw(1,:,:))./squeeze(mcal_A_rw(3,:,:)),0.84,2)';
b84_conj(index_m)=nan;
[~,~]=jbfill(ttime,a16_conj,b84_conj,rgb('darkgray'),rgb('darkgray'),0,1);
med =quantile(-squeeze(mcal_A_rw(1,:,:))./squeeze(mcal_A_rw(3,:,:)),0.5,2);
med(index_m)=nan;
plot(ttime,med,'color',rgb('darkblue'),'linewidth',1)
hold on
title('Output Growth')
set(gca,'Fontsize',gcafontsize)
hline(0,':k')
box off
grid on
ylim([0 1])
xlim([1969.75,end_date])

subplot(2,2,2)
hold on
plot(ttime,ttime*0,':r')
hold on
a16_conj = quantile(-squeeze(mcal_A_rw(2,:,:))./squeeze(mcal_A_rw(3,:,:)),0.16,2)';

a16_conj(index_m)=nan;
b84_conj = quantile(-squeeze(mcal_A_rw(2,:,:))./squeeze(mcal_A_rw(3,:,:)),0.84,2)';
b84_conj(index_m)=nan;
[~,~]=jbfill(ttime,a16_conj,b84_conj,rgb('darkgray'),rgb('darkgray'),0,1);
med =quantile(-squeeze(mcal_A_rw(2,:,:))./squeeze(mcal_A_rw(3,:,:)),0.5,2);
med(index_m)=nan;
plot(ttime,med,'color',rgb('darkblue'),'linewidth',1)
hold on
title('Core Inflation')
set(gca,'Fontsize',gcafontsize)
hline(0,':k')
box off
grid on
ylim([0 1.5])
xlim([1969.75,end_date])

subplot(2,2,3)
hold on
plot(ttime,ttime*0,':r')
hold on
a16_conj = quantile(-squeeze(mcal_A_rw(4,:,:))./squeeze(mcal_A_rw(3,:,:)),0.16,2)';
a16_conj(index_m)=nan;
b84_conj = quantile(-squeeze(mcal_A_rw(4,:,:))./squeeze(mcal_A_rw(3,:,:)),0.84,2)';
b84_conj(index_m)=nan;
[~,~]=jbfill(ttime,a16_conj,b84_conj,rgb('darkgray'),rgb('darkgray'),0,1);
med =quantile(-squeeze(mcal_A_rw(4,:,:))./squeeze(mcal_A_rw(3,:,:)),0.5,2);
med(index_m)=nan;
plot(ttime,med,'color',rgb('darkblue'),'linewidth',1)
hold on
title('M2 Growth')
set(gca,'Fontsize',gcafontsize)
hline(0,':k')
box off
grid on
ylim([0 1])
xlim([1969.75,end_date])

subplot(2,2,4)
hold on
plot(ttime,ttime*0,':r')
hold on
a16_conj = quantile(-squeeze(mcal_A_rw(5,:,:))./squeeze(mcal_A_rw(3,:,:)),0.16,2)';
a16_conj(index_m)=nan;
b84_conj = quantile(-squeeze(mcal_A_rw(5,:,:))./squeeze(mcal_A_rw(3,:,:)),0.84,2)';
b84_conj(index_m)=nan;
[~,~]=jbfill(ttime,a16_conj,b84_conj,rgb('darkgray'),rgb('darkgray'),0,1);
med =quantile(-squeeze(mcal_A_rw(5,:,:))./squeeze(mcal_A_rw(3,:,:)),0.5,2);
med(index_m)=nan;
plot(ttime,med,'color',rgb('darkblue'),'linewidth',1)
hold on
title('Credit Spread')
set(gca,'Fontsize',gcafontsize)
hline(0,':k')
box off
grid on
ylim([-4 0])
xlim([1969.75,end_date])

print([savepath, filesep, 'figure_1.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_1.png'],'-dpng');



%% Figure (2) Monetary Policy Shock
gcafontsize=14;
hFig=figure('name','Monetary policy shock');
set(hFig, 'Position', [20 20 500 350])
a16_conj = quantile(sqrt((1./squeeze(mcal_A_rw(3,:,:))).^2),0.16,2)';
a16_conj(index_m)=nan;
b84_conj = quantile(sqrt((1./squeeze(mcal_A_rw(3,:,:))).^2),0.84,2)';
b84_conj(index_m)=nan;
[~,~]=jbfill(ttime,a16_conj,b84_conj,rgb('darkgray'),rgb('darkgray'),0,0.75);
hold on
med =quantile(sqrt((1./squeeze(mcal_A_rw(3,:,:))).^2),0.5,2)
med(index_m)=nan;
plot(ttime,med,'color',rgb('darkblue'),'linewidth',1)
hold on
set(gca,'Fontsize',gcafontsize)
hline(0,':k')
box off
grid on
axis('tight')
xlim([ttime(1) ttime(end)])
ylim([0,0.8])

print([savepath, filesep, 'figure_2.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_2.png'],'-dpng');



%% Figure (3)
dates2plot=2022.25:0.25:2023.25;
hFig = figure('name','historical decomposition');
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
p1=plot(mean(squeeze(fcast_2022Q1(3,:,:)),2)*4,':k','LineWidth',2);
hold on
hold on
p2=plot(y(index_t_origin_2022Q1+1:end,3)*4,'k','LineWidth',2);
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
hold on
spf_Q1_2022.output_growth = diff(log([19893,20097,20244,20390,20533]))*400;
spf_Q1_2022.tbill3m = [0.55,0.80,1.08,1.21];
spf_Q1_2022.core_inflation = [3.1,2.5,2.3,2.3];

spf_Q2_2022.output_growth = diff(log([19736;19850;19972;20088;20193;20307]))*400;
spf_Q2_2022.tbill3m = [1;1.64;2.06;2.37;2.58];
spf_Q2_2022.core_inflation = [4.3;3.8;3.2;2.8;2.7];

p4=plot(spf_Q1_2022.tbill3m,'linestyle','none','marker','d','color','c','markerfacecolor','c');
hold on
p5=plot(spf_Q2_2022.tbill3m,'linestyle','none','marker','d','color','b','markerfacecolor','b');
leg=legend([p1 p2 h(1) h(2),p4,p5],'RC-SVAR Forecast','Data','Monetary Policy Shock','Non-Monetary Policy Shocks','SPF Q1 2022','SPF Q2 2022', 'NumColumns', 1);
set(leg,'Location','best',...
    'Orientation','vertical',...
    'NumColumns',1,...
    'FontSize',10);
set(gca,'Fontsize',12)
xtickangle(90)
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
p1=plot(mean(squeeze(fcast_2022Q1(1,:,:)),2)*4,':k','LineWidth',2);
hold on
hold on
p2=plot(y(index_t_origin_2022Q1+1:end,1)*4,'k','LineWidth',2);
hold on
p4=plot(spf_Q1_2022.output_growth,'linestyle','none','marker','d','color','c','markerfacecolor','c');
hold on
p5=plot(spf_Q2_2022.output_growth,'linestyle','none','marker','d','color','b','markerfacecolor','b');
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-5 -2.5 0 2.5 5])
ylim([-5.5,5])
ylabel('% (log-difference annualized)')
title('Output Growth')

subplot(1,3,3)
plot(mean(squeeze(fcast_2022Q1(2,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
plot(y(index_t_origin_2022Q1+1:end,2)*4,'k','LineWidth',2)
tmp=[4*mean(squeeze(shock_cont_2022Q1(2,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2022Q1(2,2:5,:,:)),1)),2)];
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
p4=plot(spf_Q1_2022.core_inflation,'linestyle','none','marker','d','color','c','markerfacecolor','c');
hold on
p5=plot(spf_Q2_2022.core_inflation,'linestyle','none','marker','d','color','b','markerfacecolor','b');
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
title('Core Inflation')
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
ylabel('% (log-difference annualized)')




print([savepath, filesep, 'figure_3.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_3.png'],'-dpng');


%% Table 1
format bank
table1=[mean(squeeze(fcast_2022Q1(3,:,:)),2)'*4;4*mean(squeeze(sum(squeeze(shock_cont_2022Q1(3,2:5,:,:)),1)),2)';4*mean(squeeze(shock_cont_2022Q1(3,1,:,:)),2)';y(index_t_origin_2022Q1+1:end,3)'*4];
display(table1)
format short



%% Figure 4
gcafontsize=12;
hFig=figure('name','counterfactuals');
set(hFig, 'Position', [0 20 950 300])

subplot(1,3,1)
a16_conj=cq16_2ppsip(3,:);
b84_conj=cq84_2ppsip(3,:);
[p1,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
a16_conj=cq16_0dot5ppsip(3,:);
b84_conj=cq84_0dot5ppsip(3,:);
[p2,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('khaki'),rgb('khaki'),0,0.55);
hold on
p3=plot(ttime(index_t_origin_2022Q1+1:end),y(end-4:end,3)*4,'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip(3,:),'linewidth',2,'color',rgb('green'))
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip(3,:),'linewidth',2,'color',rgb('yellow'))
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
a16_conj=cq16_2ppsip_ycum(1,:);
b84_conj=cq84_2ppsip_ycum(1,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
a16_conj=cq16_0dot5ppsip_ycum(1,:);
b84_conj=cq84_0dot5ppsip_ycum(1,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('khaki'),rgb('khaki'),0,0.5);
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cumsum(y(end-4:end,1)*1),'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip_ycum(1,:),'linewidth',2,'color',rgb('green'))
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip_ycum(1,:),'linewidth',2,'color',rgb('yellow'))
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
a16_conj=cq16_2ppsip(2,:);
b84_conj=cq84_2ppsip(2,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.5);
hold on
a16_conj=cq16_0dot5ppsip(2,:);
b84_conj=cq84_0dot5ppsip(2,:);
[~,~]=jbfill(ttime(index_t_origin_2022Q1+1:end),a16_conj,b84_conj,rgb('khaki'),rgb('khaki'),0,0.55);
hold on
plot(ttime(index_t_origin_2022Q1+1:end),y(end-4:end,2)*4,'linestyle','none','color',rgb('black'),'marker','o','markersize',8,'MarkerFaceColor',rgb('black'));
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_2ppsip(2,:),'linewidth',2,'color',rgb('green'));
hold on
plot(ttime(index_t_origin_2022Q1+1:end),cq50_0dot5ppsip(2,:),'linewidth',2,'color',rgb('yellow'));
axis('tight')
hold on
hline(0,':r')
box off
ylabel('%(log-difference annualized)')
grid on
title('Core Inflation')
xticks([2022.25 2022.5 2022.75 2023 2023.25])
yticks([0 2 4 6 8 10])
ylim([0,10])
xlim([2022.25,2023.25])
xticklabels({'22Q2','22Q3','22Q4','23Q1','23Q2'})
set(gca,'Fontsize',gcafontsize)
xtickangle(0)



print([savepath, filesep, 'figure_4.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_4.png'],'-dpng');




%% figure 5
[~]=pull_figure_5();


%% figure 6
dates2plot=2021.25:0.25:2021.75;
hFig = figure('name','behind the curve');
set(hFig, 'Position', [0 250 900 250])

subplot(1,3,1)
tmp=[4*mean(squeeze(shock_cont_2021Q1(3,1,1:3,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2021Q1(3,2:5,1:3,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
p1=plot(mean(squeeze(fcast_2021Q1(3,1:3,:)),2)*4,':k','LineWidth',2);
hold on
hold on
p2=plot(y(index_t_origin_2021Q1+1:index_t_origin_2021Q1+3,3)*4,'k','LineWidth',2);
xticks([find(dates2plot==2021.25) find(dates2plot==2021.5) find(dates2plot==2021.75)])
xticklabels({'2021Q2','2021Q3','2021Q4'});
box off
grid on
hold on
leg=legend([p1 p2 h(1) h(2)],'RC-SVAR Forecast','Data','Monetary Policy Shock','Non-Monetary Policy Shocks', 'NumColumns', 1);
set(leg,'Location','best',...
    'Orientation','vertical',...
    'NumColumns',1,...
    'FontSize',10);
set(gca,'Fontsize',12)
xtickangle(90)
ylabel('% (annualized)')
title('Federal Funds Rate')
legend boxoff
set(gcf, 'PaperPositionMode', 'auto');
subplot(1,3,2)
tmp=[4*mean(squeeze(shock_cont_2021Q1(1,1,1:3,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2021Q1(1,2:5,1:3,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
p1=plot(mean(squeeze(fcast_2021Q1(1,1:3,:)),2)*4,':k','LineWidth',2);
hold on
hold on
p2=plot(y(index_t_origin_2021Q1+1:index_t_origin_2021Q1+3,1)*4,'k','LineWidth',2);
hold on
xticks([find(dates2plot==2021.25) find(dates2plot==2021.5) find(dates2plot==2021.75)])
xticklabels({'2021Q2','2021Q3','2021Q4'})
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
ylabel('% (log-difference annualized)')
title('Output Growth')
subplot(1,3,3)
plot(mean(squeeze(fcast_2021Q1(2,1:3,:)),2)*4,':k','LineWidth',2)
hold on
hold on
plot(y(index_t_origin_2021Q1+1:index_t_origin_2021Q1+3,2)*4,'k','LineWidth',2)
tmp=[4*mean(squeeze(shock_cont_2021Q1(2,1,1:3,:)),2),4*mean(squeeze(sum(squeeze(shock_cont_2021Q1(2,2:5,1:3,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
xticks([find(dates2plot==2021.25) find(dates2plot==2021.5) find(dates2plot==2021.75)])
xticklabels({'2021Q2','2021Q3','2021Q4'})
title('Core Inflation')
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
ylabel('% (log-difference annualized)')


% print('results/figure_7.eps','-depsc');
% print('results/figure_7.png','-dpng');


print([savepath, filesep, 'figure_6.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_6.png'],'-dpng');

