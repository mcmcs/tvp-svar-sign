function out = plot_data(data)
  
addpath('helpfunctions/Plot')

dates = data(:,1);
Y     = data(:,2)*4; % Percent Annualized
C     = data(:,3)*4; % Percent Annualized
I     = data(:,4)*4; % Percent Annualized
P     = data(:,5)*4; % Percent Annualized
R     = data(:,6)*4; % Percent Annualized
L     = data(:,7)*1; % Percent Deviations from SS
W     = data(:,8)*4; % Percent Annualized

ftsizeaxis   = 11;
ftsizetitle  = 11;
ftsizeylabel = 11;

close all;
data_pdf = figure(1);
set(data_pdf, 'Position', [20 20 600 800])

subplot(4,2,1)
plot(dates,Y)
hline(0,'-r')
hold on
set(gca,'YTick',[-8 -4 0 4 8])
axis([dates(1) dates(end) -8 8])
ylabel('Log-diff % annualized','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('Real GDP Growth','FontSize',ftsizetitle)
 
subplot(4,2,2)
plot(dates,C)
hline(0,'-r')
hold on
set(gca,'YTick',[-4 0 4 8])
axis([dates(1) dates(end) -4 8])
ylabel('Log-diff % annualized','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('Real Consumption Growth','FontSize',ftsizetitle)

subplot(4,2,3)
plot(dates,I)
hline(0,'-r')
hold on
set(gca,'YTick',[-16 -8 0 8 16])
axis([dates(1) dates(end) -16 16])
ylabel('Log-diff % annualized','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('Real Investment Growth','FontSize',ftsizetitle)

subplot(4,2,4)
plot(dates,P)
hline(0,'-r')
hold on
set(gca,'YTick',[0 1 2 3 4 5])
axis([dates(1) dates(end) 0 5])
ylabel('Log-diff % annualized','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('GDP Deflator Inflation','FontSize',ftsizetitle)

subplot(4,2,5)
plot(dates,R)
hline(0,'-r')
hold on
set(gca,'YTick',[0 4 8 12])
axis([dates(1) dates(end) 0 12])
ylabel('% annualized','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('Nominal Interest Rate','FontSize',ftsizetitle)

subplot(4,2,6)
plot(dates,L)
hline(0,'-r')
hold on
set(gca,'YTick',[-5 -2.5 0 2.5 5])
axis([dates(1) dates(end) -5 5])
ylabel('% annualized (dev from ss)','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('Hours Worked','FontSize',ftsizetitle)

subplot(4,2,7)
plot(dates,W)
hline(0,'-r')
hold on
set(gca,'YTick',[-8 -4 0 4 8 12])
axis([dates(1) dates(end) -8 12])
ylabel('% annualized','FontSize',ftsizeylabel)
set(gca,'FontSize',ftsizeaxis)
title('Real Wage Growth','FontSize',ftsizetitle)

set(gcf, 'PaperPositionMode', 'auto');
print('figures/data','-dpdf','-bestfit')

out = data_pdf;

end

