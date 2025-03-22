clearvars;
clc;


load('flat_PC.mat');
load('steep_PC.mat');
ftsizeaxis=14;
ftsizexlabel=12;
ftsizetitle=14;
ftlinewidth = 1.0;


close all
hFig = figure('name','flat versus steep PC');
set(hFig, 'Position', [0 250 750 300])
subplot(1,2,1)
plot(oo_flat.irfs.ytildeobs_vareps_v,'color',rgb('yellow'),'LineWidth',2);
title('Output Gap')
hold on
plot(oo_steep.irfs.ytildeobs_vareps_v,'color',rgb('darkgreen'),'LineWidth',2);
title('Output Gap')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)

subplot(1,2,2)
plot(oo_flat.irfs.ppiobs_vareps_v,'color',rgb('yellow'),'LineWidth',2);
title('Inflation')
hold on
plot(oo_steep.irfs.ppiobs_vareps_v,'color',rgb('darkgreen'),'LineWidth',2);
title('Inflation')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)

% subplot(1,3,3)
% plot(oo_flat.irfs.Robs_vareps_v,'-.k');
% hold on
% plot(oo_steep.irfs.Robs_vareps_v,':r');
% legend('Flat PC','Steep PC')
% legend box off




hFig = figure('name','flat versus steep PC in levels');
set(hFig, 'Position', [0 250 750 300])
subplot(1,2,1)
plot(cumsum(oo_flat.irfs.ytildeobs_vareps_v),'color',rgb('yellow'),'LineWidth',2);
title('Output')
hold on
plot(cumsum(oo_steep.irfs.ytildeobs_vareps_v),'color',rgb('darkgreen'),'LineWidth',2);
title('Output')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)

subplot(1,2,2)
plot(cumsum(oo_flat.irfs.ppiobs_vareps_v/4),'color',rgb('yellow'),'LineWidth',2);
title('Price Level')
hold on
plot(cumsum(oo_steep.irfs.ppiobs_vareps_v/4),'color',rgb('darkgreen'),'LineWidth',2);
title('Price Level')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)