function data = get_data(start_date,end_date,sample)

currdir = pwd;


switch sample
    
    case 'GM'
        
        dates = (start_date:1/4:end_date)';
        
       
        % 1.1 GDPC1
        GDPC1    = csvread([currdir,'/data/csvfiles/GDPC1.csv'],158,1,[158 1 256 1]);
        % 1.2 GDP
        GDP      = csvread([currdir,'/data/csvfiles/GDP.csv'],158,1,[158 1 256 1]);
        % 1.3 PCND
        PCND     = csvread([currdir,'/data/csvfiles/PCND.csv'],158,1,[158 1 256 1]);
        % 1.4 PCESV
        PCESV    = csvread([currdir,'/data/csvfiles/PCESV.csv'],158,1,[158 1 256 1]);
        % 1.5 PRFI
        PRFI     = csvread([currdir,'/data/csvfiles/PRFI.csv'],158,1,[158 1 256 1]);
        % 1.6 PNFI
        PNFI     = csvread([currdir,'/data/csvfiles/PNFI.csv'],158,1,[158 1 256 1]);
        % 1.7 FEDFUNDS
        FEDFUNDS = csvread([currdir,'/data/csvfiles/FEDFUNDS.csv'],128,1,[128 1 226 1]);
        % 1.8 COMPNFB
        COMPNFB  = csvread([currdir,'/data/csvfiles/COMPNFB.csv'],158,1,[158 1 256 1]);
        % 1.9 CNP16OV
        CNP16OV  = csvread([currdir,'/data/csvfiles/CNP16OV.csv'],154,1,[154 1 252 1]);
        % 1.10 PRS85006023
        PRS85006023 = csvread([currdir,'/data/csvfiles/PRS85006023.csv'],158,1,[158 1 256 1]);
        % 1.11 CE16OV
        CE16OV = csvread([currdir,'/data/csvfiles/CE16OV.csv'],154,1,[154 1 252 1]);
        % 1.12 GDP Deflator (1.2)/(1.1)
        GDPDEF = GDP./GDPC1;
        % 1.13 Real per capita GDP = (1.1)/(1.9) 1E+9/1E+3
        GDPC1perCapita = (GDPC1./CNP16OV)*(1e9/1e3);
        % 1.14 Real per capita Consumption = (1.3 + 1.4)/((1.12)*(1.9) 1E+9/1E+3
        CperCapita = ((PCND+PCESV)./(GDPDEF.*CNP16OV))*(1e9/1e3);
        % 1.15 Real per capita Investment  = (1.5 + 1.6)/((1.12)*(1.9) 1E+9/1E+3
        XperCapita = ((PRFI+PNFI)./(GDPDEF.*CNP16OV))*(1e9/1e3);
        % 1.16 Real Wages  = (1.8)/(1.12)
        RealW = (COMPNFB./GDPDEF);
        % 1.17 Quarterly federal funds rate = 1.7/4
        FFR = FEDFUNDS/4;
        % 1.18 Labor: log(((1.10*1.11)1E3/100)/(1.9*1e3))
     
        LABOR = log(((PRS85006023.*CE16OV*1e3)/100)./(CNP16OV.*1e3))*100-mean(log(((PRS85006023.*CE16OV*1e3)/100)./(CNP16OV*1e3))*100);
        
        %% VARIABLE 1: LOG-DIFFERENCE REAL GDP PER CAPITA
        Y = diff(log(GDPC1perCapita))*100;

        %% VARIABLE 2: LOG-DIFFERENCE REAL CONSUMPTION PER CAPITA
        C = diff(log(CperCapita))*100;
        
        %% VARIABLE 3: LOG-DIFFERENCE REAL INVESTMENT PER CAPITA
        I = diff(log(XperCapita))*100;
        
        %% VARIABLE 4: LOG-DIFFERENCE GDP DEFLATOR
        P = diff(log(GDPDEF))*100;
        
        %% VARIABLE 5: FEDERAL FUNDS RATE
        R = FFR(2:end,1);
        
        %% VARIABLE 6: LABOR
        L = LABOR(2:end,1);
        
        %% VARIABLE 7: Real Wages
        W = diff(log(RealW))*100;
        
 
        data =  [dates,Y,C,I,P,R,L,W];
        
        
        
    otherwise
        
        please('Specify desired sample')
        
end

