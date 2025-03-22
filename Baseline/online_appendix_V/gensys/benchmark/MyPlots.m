clear all;
specname = 'gensys'

mkdir(['results', filesep, specname]);
savepath = [pwd, filesep, 'results', filesep, specname];



load('dynare_flat_PC.mat','oo_flat');
load('dynare_steep_PC.mat','oo_steep');
load('gensys_flat_PC.mat','gensys_flat');
load('gensys_steep_PC.mat','gensys_steep');
ftsizeaxis=14;
ftsizexlabel=12;
ftsizetitle=14;
ftlinewidth = 1.0;




hFig = figure('name','Gensys flat versus steep PC in levels');
set(hFig, 'Position', [0 250 750 300])
subplot(1,2,1)
plot(100*squeeze(cumsum(gensys_flat.irf(1,1,:))),'color',rgb('yellow'),'LineWidth',2);
title('Output')
hold on
plot(100*squeeze(cumsum(gensys_steep.irf(1,1,:))),'color',rgb('darkgreen'),'LineWidth',2);
title('Output')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)

subplot(1,2,2)
plot(100*squeeze(cumsum(gensys_flat.irf(2,1,:))),'color',rgb('yellow'),'LineWidth',2);
title('Price Level')
hold on
plot(100*squeeze(cumsum(gensys_steep.irf(2,1,:))),'color',rgb('darkgreen'),'LineWidth',2);
title('Price Level')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)


print([savepath, filesep, 'figure_PC.eps'],'-depsc');
pause(0.2);
print([savepath, filesep, 'figure_PC.png'],'-dpng');







hFig = figure('name','Compare Dynare with Gensys flat versus steep PC in levels');
set(hFig, 'Position', [0 250 750 300])
subplot(1,2,1)
plot(100*cumsum(oo_flat.irfs.ytilde_vareps_v),'color',rgb('yellow'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_flat.irf(1,1,:))),'color',rgb('yellow'),'LineWidth',2,'LineStyle',':');
title('Output')
hold on
plot(100*cumsum(oo_steep.irfs.ytilde_vareps_v),'color',rgb('darkgreen'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_steep.irf(1,1,:))),'color',rgb('darkgreen'),'LineWidth',2,'LineStyle',':');
title('Output')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)

subplot(1,2,2)
plot(100*cumsum(oo_flat.irfs.ppi_vareps_v),'color',rgb('yellow'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_flat.irf(2,1,:))),'color',rgb('yellow'),'LineWidth',2,'LineStyle',':');
title('Price Level')
hold on
plot(100*cumsum(oo_steep.irfs.ppi_vareps_v),'color',rgb('darkgreen'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_steep.irf(2,1,:))),'color',rgb('darkgreen'),'LineWidth',2,'LineStyle',':');
title('Price Level')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)





hFig = figure('name','Gensys flat versus steep PC in levels');
set(hFig, 'Position', [0 250 750 300])
subplot(1,2,1)
plot(100*cumsum(oo_flat.irfs.ytilde_vareps_v),'color',rgb('yellow'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_flat.irf(1,1,:))),'color',rgb('yellow'),'LineWidth',2,'LineStyle',':');
title('Output')
hold on
plot(100*cumsum(oo_steep.irfs.ytilde_vareps_v),'color',rgb('darkgreen'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_steep.irf(1,1,:))),'color',rgb('darkgreen'),'LineWidth',2,'LineStyle',':');
title('Output')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)

subplot(1,2,2)
plot(100*cumsum(oo_flat.irfs.ppi_vareps_v),'color',rgb('yellow'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_flat.irf(2,1,:))),'color',rgb('yellow'),'LineWidth',2,'LineStyle',':');
title('Price Level')
hold on
plot(100*cumsum(oo_steep.irfs.ppi_vareps_v),'color',rgb('darkgreen'),'LineWidth',2);
hold on
plot(100*squeeze(cumsum(gensys_steep.irf(2,1,:))),'color',rgb('darkgreen'),'LineWidth',2,'LineStyle',':');
title('Price Level')
leg=legend('Flat PC','Steep PC')
set(leg,'fontsize',12)
legend box off
box off
grid on
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)


