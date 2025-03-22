%--------------------------------------------------------------------------
% % Housekeeping
%--------------------------------------------------------------------------
clear all;
clc;
close all;



%--------------------------------------------------------------------------
% % Interpolate real GDP
%--------------------------------------------------------------------------


end_date = '2023-04-01';

% read real GDP

raw_data_GDPC1 = readtable('csvfiles/GDPC1.csv');

GDPC1          = table2array(raw_data_GDPC1(find(raw_data_GDPC1.DATE=='1959-01-01'):find(raw_data_GDPC1.DATE==end_date),2));



% read CPIAUCSL
raw_data_CPIAUCSL = readtable('csvfiles/CPIAUCSL.csv');
CPIAUCSL          = table2array(raw_data_CPIAUCSL(find(raw_data_CPIAUCSL.DATE=='1959-01-01'):find(raw_data_CPIAUCSL.DATE==end_date),2));

% read PCEPILFE
raw_data_PCEPILFE = readtable('csvfiles/PCEPILFE.csv');
PCEPILFE          = table2array(raw_data_PCEPILFE(find(raw_data_PCEPILFE.DATE=='1959-01-01'):find(raw_data_PCEPILFE.DATE==end_date),2));

% read PCEPI
raw_data_PCEPI = readtable('csvfiles/PCEPI.csv');
PCEPI          = table2array(raw_data_PCEPI(find(raw_data_PCEPI.DATE=='1959-01-01'):find(raw_data_PCEPI.DATE==end_date),2));


%--------------------------------------------------------------------------
% % Unemployment rate
%--------------------------------------------------------------------------
raw_data_UNRATE= readtable('csvfiles/UNRATE.csv');
UNRATE          = table2array(raw_data_UNRATE(find(raw_data_UNRATE.DATE=='1959-01-01'):find(raw_data_UNRATE.DATE==end_date),2));


%--------------------------------------------------------------------------
% % Federal Funds Rate
%--------------------------------------------------------------------------
raw_data_FEDFUNDS = readtable('csvfiles/FEDFUNDS.csv');
FEDFUNDS          = table2array(raw_data_FEDFUNDS(find(raw_data_FEDFUNDS.DATE=='1959-01-01'):find(raw_data_FEDFUNDS.DATE==end_date),2));



%--------------------------------------------------------------------------
% % Total reserves: Adjusted for Changes in Reserve Requirements
%--------------------------------------------------------------------------
raw_data_M1 = readtable('csvfiles/M2SL.csv');
M2          = table2array(raw_data_M1(find(raw_data_M1.DATE=='1959-01-01'):find(raw_data_M1.DATE==end_date),2));


%--------------------------------------------------------------------------
% % Corporate Credit Spread
%--------------------------------------------------------------------------
raw_data_BAA10YM = readtable('csvfiles/BAA10YM.csv');
BAA10YM          = table2array(raw_data_BAA10YM(find(raw_data_BAA10YM.DATE=='1959-01-01'):find(raw_data_BAA10YM.DATE==end_date),2));


%--------------------------------------------------------------------------
% % Export data
%--------------------------------------------------------------------------

dates = table2array(raw_data_FEDFUNDS(find(raw_data_FEDFUNDS.DATE=='1959-01-01'):find(raw_data_FEDFUNDS.DATE==end_date),1));

dataset = table(dates,GDPC1,CPIAUCSL,FEDFUNDS,M2,BAA10YM,PCEPILFE,PCEPI,UNRATE);

writetable(dataset,'csvfiles/dataset.csv')





