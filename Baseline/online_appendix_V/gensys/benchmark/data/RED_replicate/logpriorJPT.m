function [logpriorJPT,flag_ok]=logpriorJPT (param) 
% Prior for Baseline model in JPT RED 2010 
% Input PARAM is full Parameter vector, including calibrated values 
% Densities parametrized by mean and standard deviations 
% FLAG_OK is indicator == 0 if problems with prior 
prior(1)  = log(normpdf(param(1),0.30,0.05));       % 1)alpha 

prior(3)  = logBetapdf(param(3),0.5,0.15);          % 3) iotap
prior(4)  = logBetapdf(param(4),0.5,0.15);          % 4) iotaw
prior(5)  = log(normpdf(param(5),0.4,0.025));       % 5) gammastar100
prior(6)  = log(normpdf(param(6),0.2,0.025));       % 6) gammamiu100 
prior(7)  = logBetapdf( param(7),0.5,0.1);          % 7) h
prior(8)  = log(normpdf(param(8),0.15,0.05));       % 8) lambdapss
prior(9)  = log(normpdf(param(9),0.15,0.05));       % 9) lambdawss
prior(10)  = log(normpdf(param(10),0,0.5));         % 10) Lss
prior(11) = log(normpdf(param(11),0.5,0.1));        % 11) pss
prior(12) = logGammapdf(param(12),0.25,0.1);        % 12) Fbeta

prior(14) = logGammapdf(param(14),2,0.75);          % 14) niu
prior(15) = logBetapdf(param(15),0.66,0.1);         % 15) xip
prior(16) = logBetapdf(param(16),0.66,0.1);         % 16) xiw
prior(17) = logGammapdf(param(17),5,1);             % 17) chi 
prior(18) = logGammapdf(param(18),4,1);             % 18) S
prior(19) = log(normpdf(param(19),1.7,0.3));        % 19) fp
prior(20) = log(normpdf(param(20),0.125,0.05));     % 20) fy
prior(21) = logBetapdf(param(21),0.6,0.2);          % 21) rhoR
prior(22) = logBetapdf(param(22),0.4,0.2);          % 22) rhoz
prior(23) = logBetapdf(param(23),0.6,0.2);          % 23) rhog    
prior(24) = logBetapdf(param(24),0.2,0.1);          % 24) rho ISTS 
prior(25) = logBetapdf(param(25),0.6,0.2);          % 25) rholambdap
prior(26) = logBetapdf(param(26),0.6,0.2);          % 26) rholambdaw
prior(27) = logBetapdf(param(27),0.6,0.2);          % 27) rhob
prior(28) = logBetapdf(param(28),0.5,0.2);          % 28) rhoARMAlambdap
prior(29) = logBetapdf(param(29),0.5,0.2);          % 29) rhoARMAlambdaw
prior(30) = log(normpdf(param(30),0.125,0.05));     % 30) fdy

prior(32) = logBetapdf(param(32),0.6,0.2);          % 32) rho MEI 

prior(34) = log(normpdf(param(34),0.3,0.025));      % 34) GAMMA* 2
prior(35) = log(normpdf(param(35),0.6,0.025));      % 35) GAMMAMIU 2



prior(39) = logIG1pdf(param(39),0.1,1);             % sdr             
prior(40) = logIG1pdf(param(40),0.5,1);             % sdz
prior(41) = logIG1pdf(param(41),0.5,1);             % sdg
prior(42) = logIG1pdf(param(42),0.5,1);             % ISTS 
prior(43) = logIG1pdf(param(43),0.1,1);             % sdlambdap
prior(44) = logIG1pdf(param(44),0.1,1);             % sdlambdaw
prior(45) = logIG1pdf(param(45),0.1,1);             % sdb
prior(46) = logIG1pdf(param(46),0.5,1);             % MEI

if any(isnan(prior))==1 || any(isinf(prior))==1
    flag_ok=0; 
    logpriorJPT=[];
else
    flag_ok=1; 
    logpriorJPT=sum(prior);
end