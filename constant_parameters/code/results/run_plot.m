%%***************************************************************
%                         HOUSEKEEPING                          *
%%***************************************************************
clear variables;restoredefaultpath;
close all;
clc;
currdir=pwd;
cd ..
cd ..
get_help_dir_currdir=pwd;
addpath([get_help_dir_currdir,'/helpfunctions']); % set path to helper functions
cd(currdir)

load 'results.mat';

Ltilde=CLtilde;

Horizon=20;

Ltildeq50=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 50th
Ltildeq16=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 16th
Ltildeq84=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 84th

Ltildeq025=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 025th
Ltildeq975=zeros(size(Ltilde,1),size(Ltilde,2),size(Ltilde,3)); % store IRF quantile 975th


for ii=1:size(Ltilde,1)
    for jj=1:size(Ltilde,2)
        for kk=1:size(Ltilde,3)
        Ltildeq50(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.5);
        Ltildeq16(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.16);
        Ltildeq84(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.84);
        
        Ltildeq025(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.025);
        Ltildeq975(ii,jj,kk) = quantile(Ltilde(ii,jj,kk,:),0.975);
        end
    end
end

ftsizeaxis=12;
ftsizexlabel=11;
ftsizetitle=14;
ftlinewidth = 1.0;
medianwidth=2.0;

H=Horizon;



close all

hFig = figure(1);
set(hFig, 'Position', [20 20 700 350])


bands_color='silver';

subplot(2,3,1)
a=(squeeze(Ltildeq16(1:H+1,1,1)))';
b=(squeeze(Ltildeq84(1:H+1,1,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb(bands_color),rgb(bands_color),0,0.5);
hold on
plot(0:1:H,squeeze(Ltildeq50(1:H+1,1,1)),'color',rgb('black'),'LineWidth',medianwidth)
hold on
constant_para.irf_q50_output= squeeze(Ltildeq50(1:H+1,1,1));
constant_para.irf_q16_output= squeeze(Ltildeq16(1:H+1,1,1));
constant_para.irf_q84_output= squeeze(Ltildeq84(1:H+1,1,1));
hold on
hline(0,'-r')
hold on
%hold on
%[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
%set(gca,'XTick',[0;12;24;36;48;60])
%set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth )
%set(gca,'YTick',[-0.4 -0.2 0 0.2 0.4])
%axis([0 H -0.4 0.4])
%xlim([0 H])
%axis('tight')
xlabel('Months','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Output','FontSize',ftsizetitle)
box off
grid on


subplot(2,3,2)
a=(squeeze(Ltildeq16(1:H+1,2,1)))';
b=(squeeze(Ltildeq84(1:H+1,2,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb(bands_color),rgb(bands_color),0,0.5);
hold on
plot(0:1:H,squeeze(Ltildeq50(1:H+1,2,1)),'-k','LineWidth',medianwidth)
constant_para.irf_q50_cpi= squeeze(Ltildeq50(1:H+1,2,1));
constant_para.irf_q16_cpi= squeeze(Ltildeq16(1:H+1,2,1));
constant_para.irf_q84_cpi= squeeze(Ltildeq84(1:H+1,2,1));
hold on
hline(0,'-r')
hold on
%hold on
%[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
%set(gca,'XTick',[0;12;24;36;48;60])
%set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
%set(gca,'YTick',[-0.9 -0.6 -0.3 0 0.3])
%axis([0 H -0.9 0.3])
%xlabel('Months','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
%axis('tight')
xlim([0 H])
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' CPI','FontSize',ftsizetitle)
box off
grid on

subplot(2,3,3)
a=4*(squeeze(Ltildeq16(1:H+1,3,1)))';
b=4*(squeeze(Ltildeq84(1:H+1,3,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb(bands_color),rgb(bands_color),0,0.5);
hold on
plot(0:1:H,4*squeeze(Ltildeq50(1:H+1,3,1)),'-k','LineWidth',medianwidth)
constant_para.irf_q50_ffr= 4*squeeze(Ltildeq50(1:H+1,3,1));
constant_para.irf_q16_ffr= 4*squeeze(Ltildeq16(1:H+1,3,1));
constant_para.irf_q84_ffr= 4*squeeze(Ltildeq84(1:H+1,3,1));
hold on
hline(0,'-r')
hold on
%hold on
%[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
%set(gca,'XTick',[0;12;24;36;48;60])
%set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
%set(gca,'YTick',[-0.25 0 0.25 0.5])
%axis([0 H -0.25 0.5])
%axis('tight')
%xlim([0 H])
%xlabel('Months','FontSize',ftsizexlabel)
ylabel('Percentage points annualized','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Federal Funds Rate','FontSize',ftsizetitle)
box off
grid on

subplot(2,3,4)
a=(squeeze(Ltildeq16(1:H+1,4,1)))';
b=(squeeze(Ltildeq84(1:H+1,4,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb(bands_color),rgb(bands_color),0,0.5);
hold on
plot(0:1:H,squeeze(Ltildeq50(1:H+1,4,1)),'-k','LineWidth',medianwidth)
constant_para.irf_q50_m1= squeeze(Ltildeq50(1:H+1,4,1));
constant_para.irf_q16_m1= squeeze(Ltildeq16(1:H+1,4,1));
constant_para.irf_q84_m1= squeeze(Ltildeq84(1:H+1,4,1));
hold on
hline(0,'-r')
hold on
%hold on
%[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
%set(gca,'XTick',[0;12;24;36;48;60])
%set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
%set(gca,'YTick',[-1.5 -1 -0.5 0])
%axis([0 60 -1.2 1.2])
%axis('tight')
xlim([0 H])
%xlabel('Months','FontSize',ftsizexlabel)
ylabel('Percent','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' M1','FontSize',ftsizetitle)
box off
grid on

subplot(2,3,5)
a=4*(squeeze(Ltildeq16(1:H+1,5,1)))';
b=4*(squeeze(Ltildeq84(1:H+1,5,1)))';
x = 0:1:H;
[~,~]=jbfill(x,a,b,rgb(bands_color),rgb(bands_color),0,0.5);
hold on
plot(0:1:H,4*squeeze(Ltildeq50(1:H+1,5,1)),'-k','LineWidth',medianwidth)
constant_para.irf_q50_m1= 4*squeeze(Ltildeq50(1:H+1,5,1));
constant_para.irf_q16_m1= 4*squeeze(Ltildeq16(1:H+1,5,1));
constant_para.irf_q84_m1= 4*squeeze(Ltildeq84(1:H+1,5,1));
hold on
hline(0,'-r')
hold on
%hold on
%[~,~]=jbfill(x,a95,b95,rgb('royalblue'),rgb('royalblue'),0,0.5);
%set(gca,'XTick',[0;12;24;36;48;60])
%set(gca,'XTickLabel',['0 ';' 1';' 2';' 3';' 4';' 5'])
set(gca,'LineWidth',ftlinewidth)
%set(gca,'YTick',[-1.5 -1 -0.5 0])
%axis([0 60 -1.2 1.2])
%axis('tight')
xlim([0 H])
%xlabel('Months','FontSize',ftsizexlabel)
ylabel('Percentage points annualized','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(' Corporate credit spread','FontSize',ftsizetitle)
box off
grid on


set(gcf, 'PaperPositionMode', 'auto');




print -dpng  'uhlig_2017_constant_para.png'
print('uhlig_2017_constant_para.eps','-depsc');


save('constant_para_irfs.mat','constant_para');








%% Historical decomposition
t_origin = 2022;
ttime    = 1960.25:0.25:2023.25;

%
index_t_origin=find(ttime==t_origin);



H = 4;
info_nlag=nlag;
nd =size(Ltilde,4);
info_nvar=size(Ltilde,2);
info_m = size(Aplustilde,1);
nvar=info_nvar;
fcast      = nan(info_nvar,H+1,nd);
shock_cont = nan(info_nvar,info_nvar,H+1,nd);
shocks_RC_fixed_para=nan(info_nvar,H+1,nd);
L0_check = zeros(nd,1);
MPshocks_RC  = nan(info_nvar,nd,length(ttime));
neg_forecast_index = zeros(nd,1);
e=eye(info_nvar);

y=Y;
x=X;

for i=1:1:nd



    
    

    A = A0tilde(:,:,i); %reshape(mcal_A_tmp{i}(:,index_t_origin),info_nvar,info_nvar);
    F = Aplustilde(:,:,i); %reshape(mcal_F_tmp{i}(:,index_t_origin),info_m,info_nvar);
  

    structpara = [A(:);F(:)];
    
 
    
    LIRF       = IRF_horizons(structpara, info_nvar, info_nlag, info_m, 0:5);


    B = F/A;
    
    for t=1:length(ttime)
    MPshocks_RC(:,i,t) = ((y(t,:)-x(t,:)*B)*A)';
    end
%         non_stationary = 0;
%         Fmatrix = zeros(nvar*nlag,nvar*nlag);
%         for ell=1:nlag
%             F(1:nvar,1 + (ell-1)*nvar:ell*nvar)=B(1+(ell-1)*nvar:ell*nvar,:)';
%             if ell<nlag
%                 Fmatrix(1+ell*nvar:nvar+ell*nvar,1+(ell-1)*nvar:ell*nvar)=eye(nvar);
%             end
%         end
%         EF=eig(Fmatrix);
%         for ii=1:nvar*nlag
%             if norm(EF(ii))>1
%                 non_stationary = 1;
%             end
%         end
% if  non_stationary==1
% keyboard
% end

    L0_check(i,1) = LIRF(1,1);

    for t=1:(H+1)

        % y(t)' A(t) = x(t)' F(t) + eps(t)'
        % y(t)' = x(t)' F(t) (A(t))^-1 + eps(t)' A(t)^-1
        % y(t)' = x(t)' B(t) + eps(t)' A(t)^-1
        % y(t)' = x(t)' B(t) + u(t)'
        % u(t)' = y(index_t_origin+t,:)-x(index_t_origin+t,:)*B
        % u(t)' = eps(t)' A(t)^-1
        % u(t)'*A = eps(t)' 
       
        shocks_RC_fixed_para(:,t,i) = ((y(index_t_origin+t,:)-x(index_t_origin+t,:)*B)*A)';
        
    end
    

    for i_var=1:info_nvar
        for j_shock=1:info_nvar

            for h=0:H

                switch h
                    case 0
                        xtprime = [y(index_t_origin,:),y(index_t_origin-1,:),y(index_t_origin-2,:),y(index_t_origin-3,:),1];
                    case 1
                        xtprime = [squeeze(fcast(:,h,i))',y(index_t_origin,:),y(index_t_origin-1,:),y(index_t_origin-2,:),1];
                   case 2
                        xtprime = [squeeze(fcast(:,h,i))',squeeze(fcast(:,h-1,i))',y(index_t_origin,:),y(index_t_origin-1,:),1];
                   case 3
                        xtprime = [squeeze(fcast(:,h,i))',squeeze(fcast(:,h-1,i))',squeeze(fcast(:,h-2,i))',y(index_t_origin,:),1];
                    otherwise
                        xtprime = [squeeze(fcast(:,h,i))',squeeze(fcast(:,h-1,i))',squeeze(fcast(:,h-2,i))',squeeze(fcast(:,h-3,i))',1];
                end
                fcast(:,h+1,i) = B'*xtprime';

                  if min(fcast(:,h+1,i))<0
                    neg_forecast_index(i,1)=1;
                end

                tmp=0;

                for ell=0:h
                    tmp= tmp + e(:,i_var)'*LIRF(1+ell*info_nvar:(ell+1)*info_nvar,:)*e(:,j_shock)*e(:,j_shock)'*shocks_RC_fixed_para(:,1+h-ell,i);
                end
                   shock_cont(i_var,j_shock,h+1,i) = tmp;
            end
        end

    end

end



close all

spf_Q1_2022.output_growth = diff(log([19893,20097,20244,20390,20533]))*400;
spf_Q1_2022.output_growth_2023 = diff(log([20699,21185]))*100;
spf_Q1_2022.tbill3m = [0.55,0.80,1.08,1.21];
spf_Q1_2022.core_inflation = [3.1,2.5,2.3,2.3];

spf_Q2_2022.output_growth = diff(log([19736;19850;19972;20088;20193;20307]))*400;
spf_Q2_2022.output_growth_2023 = diff(log([20356,20765]))*100;
spf_Q2_2022.tbill3m = [1;1.64;2.06;2.37;2.58];
spf_Q2_2022.core_inflation = [4.3;3.8;3.2;2.8;2.7];


close all
%% historical decomposition 1 x 3
dates2plot=2022.25:0.25:2023.25;
hFig = figure(2);
set(hFig, 'Position', [0 250 1000 250])

subplot(1,3,1)

tmp=[4*mean(squeeze(shock_cont(3,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont(3,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%h=bar(shocks(index_bd_shock_decomposition,[1:7, 11],1).*(shocks(index_bd_shock_decomposition,[1:7, 11],1)<0),'Stacked');
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%p1=plot(mean(squeeze(fcast(3,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
%p2=plot(y(index_t_origin+1:end,3)*4,'k','LineWidth',2)
%vline(find(dates2plot==2022),':k')
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
%leg=legend('Forecast as of 2021Q4','Data','Monetary Policy Shock','Sum of Other Shocks', 'NumColumns', 4);
% leg=legend('','','','', 'NumColumns', 4);
% set(leg,'Fontsize',12)
% set(leg,'location','bestoutside')
% set(leg,'orientation','horizontal')
%legend boxoff
box off
grid on
hold on
%p3=plot(sep_march_2022.ffr,'linestyle','none','marker','o','color','g','markerfacecolor','g')
%hold on
%p4=plot(spf_Q1_2022.tbill3m,'linestyle','none','marker','d','color','c','markerfacecolor','c')
hold on
%p5=plot(spf_Q2_2022.tbill3m,'linestyle','none','marker','d','color','b','markerfacecolor','b')
%leg=legend([p1 p2 h(1) h(2),p4,p5],'RC-SVAR Forecast','Data','Monetary Policy Shock','Non-Monetary Policy Shocks','SPF Q1 2022','SPF Q2 2022', 'NumColumns', 1);
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

tmp=[4*mean(squeeze(shock_cont(1,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont(1,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%h=bar(shocks(index_bd_shock_decomposition,[1:7, 11],1).*(shocks(index_bd_shock_decomposition,[1:7, 11],1)<0),'Stacked');
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%p1=plot(mean(squeeze(fcast(1,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
%p2=plot(y(index_t_origin+1:end,1)*4,'k','LineWidth',2)
hold on
%p3=plot(sep_march_2022.output_growth,'linestyle','none','marker','o','color','g','markerfacecolor','g')
hold on
%p4=plot(spf_Q1_2022.output_growth,'linestyle','none','marker','d','color','c','markerfacecolor','c')
hold on
%p5=plot(spf_Q2_2022.output_growth,'linestyle','none','marker','d','color','b','markerfacecolor','b')
%vline(find(dates2plot==2022),':k')
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
%leg=legend('','','','', 'NumColumns', 4);
%set(leg,'Fontsize',12)
%set(leg,'location','best')
%set(leg,'orientation','horizontal')
%legend boxoff
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-3 -2 -1 0 1 2 3 4])
ylim([-3,4])
ylabel('% (log-difference annualized)')
title('Output Growth')

subplot(1,3,3)
%plot(mean(squeeze(fcast(2,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
%plot(y(index_t_origin+1:end,2)*4,'k','LineWidth',2)

tmp=[4*mean(squeeze(shock_cont(2,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont(2,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%h=bar(shocks(index_bd_shock_decomposition,[1:7, 11],1).*(shocks(index_bd_shock_decomposition,[1:7, 11],1)<0),'Stacked');
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
hold on
%p3=plot(sep_march_2022.core_inflation,'linestyle','none','marker','o','color','g','markerfacecolor','g')
hold on
%p4=plot(spf_Q1_2022.core_inflation,'linestyle','none','marker','d','color','c','markerfacecolor','c')
hold on
%p5=plot(spf_Q2_2022.core_inflation,'linestyle','none','marker','d','color','b','markerfacecolor','b')
%vline(find(dates2plot==2022),':k')
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
% leg=legend('Forecast as of 2021Q4','Data','Monetary Policy Shock','Sum of Other Shocks', 'NumColumns', 4);
% set(leg,'Fontsize',12)
% set(leg,'location','best')
% set(leg,'orientation','horizontal')
%legend boxoff
title('Core Inflation')
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-1.5 -0.5 0 0.5 1 1.5 2])
ylim([-1.5,2])
ylabel('% (log-difference annualized)')



print('shock_decomp_1by3_sys_uhlig_spf_march_june_CONSTANT_PARA_R1.eps','-depsc');
print('shock_decomp_1by3_sys_uhlig_spf_march_june_CONSTANT_PARA_R1.png','-dpng');





close all
%% historical decomposition 1 x 3
dates2plot=2022.25:0.25:2023.25;
hFig = figure(2);
set(hFig, 'Position', [0 250 1000 250])

subplot(1,3,1)

tmp=[4*mean(squeeze(shock_cont(3,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont(3,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%h=bar(shocks(index_bd_shock_decomposition,[1:7, 11],1).*(shocks(index_bd_shock_decomposition,[1:7, 11],1)<0),'Stacked');
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%p1=plot(mean(squeeze(fcast(3,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
%p2=plot(y(index_t_origin+1:end,3)*4,'k','LineWidth',2)
%vline(find(dates2plot==2022),':k')
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
%leg=legend('Forecast as of 2021Q4','Data','Monetary Policy Shock','Sum of Other Shocks', 'NumColumns', 4);
% leg=legend('','','','', 'NumColumns', 4);
% set(leg,'Fontsize',12)
% set(leg,'location','bestoutside')
% set(leg,'orientation','horizontal')
%legend boxoff
box off
grid on
hold on
%p3=plot(sep_march_2022.ffr,'linestyle','none','marker','o','color','g','markerfacecolor','g')
%hold on
%p4=plot(spf_Q1_2022.tbill3m,'linestyle','none','marker','d','color','c','markerfacecolor','c')
hold on
%p5=plot(spf_Q2_2022.tbill3m,'linestyle','none','marker','d','color','b','markerfacecolor','b')
%leg=legend([p1 p2 h(1) h(2),p4,p5],'RC-SVAR Forecast','Data','Monetary Policy Shock','Non-Monetary Policy Shocks','SPF Q1 2022','SPF Q2 2022', 'NumColumns', 1);
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

tmp=[4*mean(squeeze(shock_cont(4,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont(4,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%h=bar(shocks(index_bd_shock_decomposition,[1:7, 11],1).*(shocks(index_bd_shock_decomposition,[1:7, 11],1)<0),'Stacked');
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%p1=plot(mean(squeeze(fcast(1,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
%p2=plot(y(index_t_origin+1:end,1)*4,'k','LineWidth',2)
hold on
%p3=plot(sep_march_2022.output_growth,'linestyle','none','marker','o','color','g','markerfacecolor','g')
hold on
%p4=plot(spf_Q1_2022.output_growth,'linestyle','none','marker','d','color','c','markerfacecolor','c')
hold on
%p5=plot(spf_Q2_2022.output_growth,'linestyle','none','marker','d','color','b','markerfacecolor','b')
%vline(find(dates2plot==2022),':k')
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
%leg=legend('','','','', 'NumColumns', 4);
%set(leg,'Fontsize',12)
%set(leg,'location','best')
%set(leg,'orientation','horizontal')
%legend boxoff
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-12 -9 -6 -3 0])
ylim([-12,0])
ylabel('% (log-difference annualized)')
title('Money Growth')

subplot(1,3,3)

tmp=[4*mean(squeeze(shock_cont(5,1,:,:)),2),4*mean(squeeze(sum(squeeze(shock_cont(5,2:5,:,:)),1)),2)];
h=bar(tmp.*(tmp>0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%h=bar(shocks(index_bd_shock_decomposition,[1:7, 11],1).*(shocks(index_bd_shock_decomposition,[1:7, 11],1)<0),'Stacked');
h=bar(tmp.*(tmp<0),'Stacked');
set(h(1),{'FaceColor'},{'r'});
set(h(2),{'FaceColor'},{'y'});
hold on
%p1=plot(mean(squeeze(fcast(1,:,:)),2)*4,':k','LineWidth',2)
hold on
hold on
%p2=plot(y(index_t_origin+1:end,1)*4,'k','LineWidth',2)
hold on
%p3=plot(sep_march_2022.output_growth,'linestyle','none','marker','o','color','g','markerfacecolor','g')
hold on
%p4=plot(spf_Q1_2022.output_growth,'linestyle','none','marker','d','color','c','markerfacecolor','c')
hold on
%p5=plot(spf_Q2_2022.output_growth,'linestyle','none','marker','d','color','b','markerfacecolor','b')
%vline(find(dates2plot==2022),':k')
xticks([find(dates2plot==2022.25) find(dates2plot==2022.5) find(dates2plot==2022.75) find(dates2plot==2023) find(dates2plot==2023.25)])
xticklabels({'2022Q2','2022Q3','2022Q4','2023Q1','2023Q2'})
%leg=legend('','','','', 'NumColumns', 4);
%set(leg,'Fontsize',12)
%set(leg,'location','best')
%set(leg,'orientation','horizontal')
%legend boxoff
box off
grid on
set(gca,'Fontsize',12)
xtickangle(90)
yticks([-0.4 -0.2 0 0.2 0.4])
ylim([-0.4,0.4])
ylabel('% (annualized)')
title(' Credit Spreads')

print('shock_decomp_1by3_sys_uhlig_spf_march_june_CONSTANT_PARA_R1_MoneyGrowth.eps','-depsc');
print('shock_decomp_1by3_sys_uhlig_spf_march_june_CONSTANT_PARA_R1_MoneyGrowth.png','-dpng');



% %% implied shocks
% 
% nburn=0;
% 
% close all
% gcafontsize=14;
% hFig=figure('name','shocks');
% set(hFig, 'Position', [0 20 450 300])
% ttime_start=find(ttime==2021);
% a16_conj = quantile(squeeze(MPshocks_RC(1,nburn+1:end,ttime_start:end))',0.16,2)';
% b84_conj = quantile(squeeze(MPshocks_RC(1,nburn+1:end,ttime_start:end))',0.84,2)';
% 
% a16_conj_baseline = a16_conj;
% b84_conj_baseline = b84_conj;
% 
% [~,~]=jbfill(ttime(ttime_start:end),a16_conj,b84_conj,rgb('lightgreen'),rgb('lightgreen'),0,0.75);
% hold on
% med_baseline =median(squeeze(MPshocks_RC(1,nburn+1:end,ttime_start:end)),1);
% plot(ttime(ttime_start:end),median(squeeze(MPshocks_RC(1,nburn+1:end,ttime_start:end)),1),'color',rgb('darkgreen'),'linewidth',1,'LineStyle','-')
% hold on
% hline(0,':k')
% hold on
% hold on
% %xline(ttime(ttime==2022.5),'--r','2022Q3 (Contractionary)')
% xlim([ttime(ttime_start) ttime(end)])
% hold on
%  recessionplot
%  hold on
% %title('Monetary Policy Shocks')
% box off
% grid on
% xticks([ttime(ttime_start) ttime(ttime_start+2) ttime(ttime_start+4)  ttime(ttime_start+6) ttime(ttime_start+8)])
% grid on
% xticklabels({'2021Q1','2021Q3','2022Q1','2022Q3','2023Q1'})
% set(gca,'Fontsize',gcafontsize)
% print('sys_and_uhlig_shocks_zoom_2021Q1_2023Q2_CONSTANT_PARA_R1.png','-dpng');
% print('sys_and_uhlig_shocks_zoom_2021Q1_2023Q2_CONSTANT_PARA_R1.eps','-depsc');
% 
% 
% save('figure_sys_and_uhlig_shocks_zoom_2021Q1_2023Q2','med_baseline','a16_conj_baseline','b84_conj_baseline')
% 
