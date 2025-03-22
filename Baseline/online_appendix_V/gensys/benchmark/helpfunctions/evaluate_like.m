function [loglh] = evaluate_like(ttheta,Y,T1,T0)
% This function evaluates the likelihood. The code is based on the Gauss procedure evalcgg by Lubik and Schorfheide (2004)


% 0. Preliminaries
cal_para    = get_calibrated_parameters();
Ppiss       = cal_para(8,1);
Ztildess    = cal_para(9,1);
bbeta       = cal_para(11,1);
Rss         = Ztildess*Ppiss/bbeta; 

% 1. Define matrices in Hamilton (1994)'s notation (Page 377)
% \xi_{t+1} = F \xi_{t} + v_{t}
% y_{t}     = A' x_{t} + H' s_{t} + w_{t}
% E(v_{t}v_{tau}^{\prime}) = Q if t=ttau
%                            0 Otherwise
% E(w_{t}w_{tau}^{\prime}) = R if t=ttau
%                            0 Otherwise

% Model notation
% s_{t+1} = T1 s_{t} + T0 eps_{t}
% y_{t}   = A' x_{t} + H' s_{t} + w_{t}

% Mapping Hamilton's notation to model
% s_{t+1} = \xi_{t+1}
% v_{t}   = T0 * eps_{t}
% x_{t}   = 1
% T1      = F
% Q       = T0*Oomega*T0'

nvars        = size(Y,2);
nstates      = size(T1,1);
 
At          = zeros(nvars,1);
At(1,1)     = log(Ztildess)*100;
At(2,1)     = log(Ztildess)*100;
At(3,1)     = log(Ztildess)*100;
At(4,1)     = 100*log(Ppiss);
At(5,1)     = 100*log(Rss);
At(6,1)     = 0;
At(7,1)     = log(Ztildess)*100;
A           = At';


Ht          = zeros(nvars,nstates);
Ht(1,25)    = 100; Ht(1,30) = -100; Ht(1,7) = 100;
Ht(2,2)     = 100; Ht(2,32) = -100; Ht(2,7) = 100;
Ht(3,6)     = 100; Ht(3,33) = -100; Ht(3,7) = 100;
Ht(4,20)    = 100;
Ht(5,11)    = 100;
Ht(6,3)     = 100;
Ht(7,1)     = 100; Ht(7,34) = -100; Ht(7,7) = 100;
H           = Ht';

F           = T1;


Oomega      = zeros(nvars,nvars);
Oomega(1,1) = (ttheta(18))^2;
Oomega(2,2) = (ttheta(19))^2;
Oomega(3,3) = (ttheta(20))^2;
Oomega(4,4) = (ttheta(21))^2;
Oomega(5,5) = (cal_para(13,1))^2;
Oomega(6,6) = (ttheta(22))^2;
Oomega(7,7) = (ttheta(23))^2;

Q           = T0*Oomega*T0';
Q           = 0.5*(Q+Q');


R        = zeros(nvars,nvars);
R(4,4)   = (100*ttheta(24))^2;
R(7,7)   = (100*ttheta(25))^2;




% Set loglh=-1e8 if Q is not 
if max(eig(Q))<=0
   loglh=-1e8;
   return
end


% 2. Set initial values for the kalman filter
Xxi      = zeros(nstates,1); %Xxi denotes State vector
vecQ     = Q(:);
vecP     = (eye(nstates*nstates)-kron(F,F))\vecQ;
P        = reshape(vecP,size(F,1),size(F,1));
P        = (P+P')/2;


% 3. Initialize iteration
t=1;
T     = size(Y,1);
loglh = 0;

while t<=T; 
    

y = Y(t,:)'; % y is a n by 1 vector



% Forecasting
Xxihat = F*Xxi;           % 13.2.20
Phat   = F*P*F' + Q;      % 13.2.21

Ssigmay = H'*Phat*H + R;
Ssigmay = 0.5*(Ssigmay+Ssigmay');
yhat    = A' + H'*Xxihat;
ferror  = y-yhat;


% Set loglh=-1e8 Ssigmay is not well-conditioned
if rcond(Ssigmay)<1e-8
   loglh=-1e8;
   display('Ssigmay is not well-conditioned')
   return
end

Term1 = -0.5*nvars*log(2*pi)-0.5*log(det(Ssigmay));
Term2 = -0.5*ferror'*(Ssigmay\ferror);


loglh = Term1 + Term2 + loglh;

% Updating 
Xxi = Xxihat + Phat*H*(Ssigmay\ferror);     %13.2.15
P   = Phat   - Phat*H*(Ssigmay\(H'*Phat')); %13.2.16

t   = t + 1;

end

