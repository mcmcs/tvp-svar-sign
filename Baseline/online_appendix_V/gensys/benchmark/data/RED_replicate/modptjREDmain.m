 function [Gone,Cone,RRone,euone,SDXone,zmat,NY,NX,NETA,gev,ncof,...
      Gtwo,RRtwo,eutwo,SDXtwo,Ctwo,Tbreak]...
     =modptjREDmain(estpr,calpr,indp,con,ind_sv)
% =========================================================================
% MOD PTJ RED main 
% 
% Split sample estimation of baseline mode RED 2010 
% For given parameter values, this function 
% 1) puts the model of Justiniano, Primiceri and Tambalotti (2010)
%    in Gensys' canonical form
% 2) solves the RE systemS of equations using Chris Sims' Gensys
% 3) augments the resulting state evolution equation with the observation
%    equation
% The solution for each subsample takes the form:  x(t) = G<i> * x(t-1) + R<i>* e(t)
%                               y(t) = zmat * x(t) + C<i>
%                               e(t) ~ N(0,SDX<i>'*SDX<i>)  i=1,2 
%
% Tbreak is the date of break (observation number)
%
% Alejandro Justiniano, Giorgio Primiceri and Andrea Tambalotti 01/30/09  
% =========================================================================

y=1;            % output
k=2;            % capital
L=3;            % hours
Rk=4;           % return on capital
w=5;            % real wages
p=6;            % inflation
s=7;            % marginal cost
lambda=8;       % multiplier
c=9;            % consumption
R=10;           % interst rate
u=11;           % capital utilization
phi=12;         % multiplier
i=13;           % investment
kbar=14;        % "gross" capital
wtilda=15;      % auxiliary variable for wage equation
x=16;           % marginal utility of labor
z=17;           % productivity shock
g=18;           % public spending shock
miu=19;         % relative price (non stationary)
lambdap=20;     % mark-up shock
psi=21;         % preference shock
b=22;           % discount factor shock 
ep=23;          % expectational variables
ec=24;
elambda=25;
ephi=26;
eRk=27;
ei=28;
ewtilda=29;
yobs=30; 
cobs=31;
iobs=32; 
wobs=33; 
ARMAlambdap=34;
ARMApsi=35;
mp=36;

% variables for the flex prices and wages version of the model
%===
ystar=36+1;            % output
kstar=36+2;            % capital
Lstar=36+3;            % hours
Rkstar=36+4;           % return on capital
wstar=36+5;            % real wages
sstar=36+6;            % marginal cost
lambdastar=36+7;       % multiplier
cstar=36+8;            % consumption
Rstar=36+9;           % interst rate
ustar=36+10;           % capital utilization
phistar=36+11;         % multiplier
istar=36+12;           % investment
kbarstar=36+13;        % "gross" capital
wtildastar=36+14;      % auxiliary variable for wage equation
xstar=36+15;           % marginal utility of labor
ecstar=36+16;
elambdastar=36+17;
ephistar=36+18;
eRkstar=36+19;
eistar=36+20;
ygaplev=36+21; 
dagtrend=36+22;  % First difference of the aggregate trend 

upsilon=36+23;  % new, stationary, investment shock
gdp    =36+24; 
gdpstar=36+25; 

NY=36+25; 

% =====
% INDEX for Shocks  
% ****************************************
NX=8; 

Rs=1;           % MP
zs=2;           % Technology 
gs=3;           % Public spending
mius=4;         % Investment     
lambdaps=5;     % Mark-up
psis=6;         % Leisure preference
bs=7;           % Discount factor
upsilons=8;     % investment shock


% Eta - Expectational errors
% ---------------------------
NETA=7+5; 

pex=1; 
cex=2;
lambdaex=3;
phiex=4;
Rkex=5;
iex=6;
wtildaex=7;

%=============
cexstar=8;
lambdaexstar=9;
phiexstar=10;
Rkexstar=11;
iexstar=12;
%============

% Index for the parameters 
% -----------------------
numpar=      46;  % Number of parameters 
ncof=        38;  % Number of coefficients not corresponding to 
                  % standard deviations 
spadd=5; 

if nargin < 5 || isempty( ind_sv ) 
    flag_sv = 0 ;
else 
    flag_sv = 1 ;
end 


if nargin == 1 || isempty( calpr ) 
    param=estpr  ; 
else
    param=zeros(numpar,1); 
    param(indp)=estpr; 
    param(con)=calpr; 
end; 


% Non SD parameters 
% -------------------
alpha=param(1);     % share of capital in the prod. function
delta=param(2);     % capital depreciation rate
iotap=param(3);     % price indexation (=0 is static indexation, =1 is dynamic)
iotaw=param(4);     % wages indexation
gammastar100=param(5);  % steady state growth rate of technology
gammamiu100=param(6);   % steady state growth rate of capital embodied technology
h=param(7);         % habit formation parameter
lambdapss=param(8); % steady state of mark-up shock (pins down steady state of wages)
lambdaw=param(9);   % wage mark-up
Lss=param(10);      % steady state for leisure (the ss for leisure is pinned down by psiss. convenient to parameterize the model in terms of Lss)
pss100=param(11);   % steady state quarterly inflation (multiplied by 100)
Fbeta=param(12);    % weird parameter of SW that is a function of beta and the steady state quarterly real rate of interest (ensures that beta is less than 1)
gss=param(13);      % steady state of the public spending shock
niu=param(14);      % curvature of the utility function for leisure
xip=param(15);      % price stickiness
xiw=param(16);      % wage stickiness
chi=param(17);      % elasticity of the capital utilization cost function
S=param(18);        % investment adjustment cost
fp=param(19);       % reaction to inflation
fy=param(20);       % reaction to output gap
rhoR=param(21); rhoz=param(22); rhog=param(23); rhomiu=param(24); rholambdap=param(25); rhopsi=param(26); rhob=param(27); 
rhoARMAlambdap=param(28); rhoARMApsi=param(29);
fdy=param(30);      % reaction to output gap growth (like in SW)
rhomp=param(31);    % autocorrelation of MP shock
rhoupsilon=param(32); 

flag_relpiobs=param(33); 
if flag_relpiobs~=1 && flag_relpiobs~=0; 
    error('FLAG RELPIOBS must be 0 or 1') 
end 

% ================================
% Parameters for the second sample 
% ================================
gammastar100s2=param(34); 
gammamiu100s2 =param(35); 
rhomius2      =param(36); 
if rhomius2==999; 
    rhomius2=rhomiu; 
end 
sdmius2      =param(37); 
if sdmius2==999; 
    sdmius2=param(42); 
end 
Tbreak       =param(38); 

% Parameters for the second sample 
partwo=zeros(numpar-spadd,1); 
partwo(1:4)=param(1:4); 
partwo(5)  =gammastar100s2; 
partwo(6)  =gammamiu100s2 ; 
partwo(7:33)=param(7:33); 
partwo(24  )=rhomius2; 
partwo(34:end)=param(39:end); 
partwo(37   )=sdmius2; 

% STEADY STATE
% ----------------------------
gammastar=gammastar100/100;
gammamiu=gammamiu100/100;
beta=100/(Fbeta+100);
rss=exp(gammastar)/beta-1;
rss100=rss*100;
pss=pss100/100;
gss=1/(1-gss);

Rkss=exp(gammastar+gammamiu)/beta-1+delta;
sss=1/(1+lambdapss);
wss=(sss*((1-alpha)^(1-alpha))/((alpha^(-alpha))*Rkss^alpha))^(1/(1-alpha));
kLss=(wss/Rkss)*alpha/(1-alpha);
FLss=(kLss^alpha-Rkss*kLss-wss);
yLss=kLss^alpha-FLss;
Lsscale=1;     
kss=kLss*Lsscale;
iss=(1-(1-delta)*exp(-gammastar-gammamiu))*kss*exp(gammastar+gammamiu);
F=FLss*Lsscale;
yss=yLss*Lsscale;
css=yss/gss-iss;    
    
% System Matrices 
GAM0 = zeros(NY,NY) ;
GAM1 = zeros(NY,NY) ;
PSI = zeros(NY,NX) ;
PPI = zeros(NY,NETA) ;
C = zeros(NY,1) ;


% eq 1, prod. function (y)
% ----------------------------------------
GAM0(y,y)=1;
GAM0(y,k)=-((yss+F)/yss)*alpha;
GAM0(y,L)=-((yss+F)/yss)*(1-alpha);

% ===
GAM0(ystar,ystar)=1;
GAM0(ystar,kstar)=-((yss+F)/yss)*alpha;
GAM0(ystar,Lstar)=-((yss+F)/yss)*(1-alpha);
%===

% GDP 
GAM0(gdp,gdp)= -1; 
GAM0(gdp,y)=1; 
GAM0(gdp,u)= -kss*Rkss/yss; 

% GDPSTAR 
GAM0(gdpstar,gdpstar)= -1; 
GAM0(gdpstar,ystar)=1; 
GAM0(gdpstar,ustar)= -kss*Rkss/yss;

% eq 2, cost minimization (L)
% ------------------------------------------
GAM0(L,Rk)=1;
GAM0(L,w)=-1;
GAM0(L,k)=1;
GAM0(L,L)=-1;

% ===
GAM0(Lstar,Rkstar)=1;
GAM0(Lstar,wstar)=-1;
GAM0(Lstar,kstar)=1;
GAM0(Lstar,Lstar)=-1;
% ===

% eq 4, Phillips Curve (p)
% ------------------------------------------
GAM0(p,p)=1;
GAM0(p,ep)=-beta/(1+iotap*beta);
GAM1(p,p)=iotap/(1+iotap*beta);
GAM0(p,s)=-((1-beta*xip)*(1-xip)/((1+iotap*beta)*xip));
GAM0(p,lambdap)=-1;        
%===
GAM0(Rstar,sstar)=1;
%===

% eq 5, marginal cost (s)
% --------------------------------------------
GAM0(s,s)=1;
GAM0(s,Rk)=-alpha;
GAM0(s,w)=-(1-alpha);

%===
GAM0(sstar,sstar)=1;
GAM0(sstar,Rkstar)=-alpha;
GAM0(sstar,wstar)=-(1-alpha);
%===

% eq 7, consumption foc (c)
% --------------------------------------------
expg=exp(gammastar);
GAM0(c,lambda)=(expg-h*beta)*(expg-h);
GAM0(c,b)=-(expg-h*beta*rhob)*(expg-h)/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)]; 
GAM0(c,z)=-(beta*h*expg*rhoz-h*expg);
GAM0(c,miu)=[-(beta*h*expg*rhomiu-h*expg)]*alpha/(1-alpha);
GAM0(c,c)=(expg^2+beta*h^2);
GAM0(c,ec)=-beta*h*expg;
GAM1(c,c)=expg*h;

%===
expg=exp(gammastar);
GAM0(cstar,lambdastar)=(expg-h*beta)*(expg-h);
GAM0(cstar,b)=-(expg-h*beta*rhob)*(expg-h)/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];  
GAM0(cstar,z)=-(beta*h*expg*rhoz-h*expg);
GAM0(cstar,miu)=[-(beta*h*expg*rhomiu-h*expg)]*alpha/(1-alpha);
GAM0(cstar,cstar)=(expg^2+beta*h^2);
GAM0(cstar,ecstar)=-beta*h*expg;
GAM1(cstar,cstar)=expg*h;
%===

% eq 8, euler equation (lambda)
% --------------------------------------------
GAM0(lambda,lambda)=1;
GAM0(lambda,R)=-1;
GAM0(lambda,elambda)=-1;
GAM0(lambda,ep)=1;
GAM0(lambda,z)=rhoz;
GAM0(lambda,miu)=[rhomiu]*alpha/(1-alpha);

%===
GAM0(lambdastar,lambdastar)=1;
GAM0(lambdastar,Rstar)=-1;
GAM0(lambdastar,elambdastar)=-1;
GAM0(lambdastar,z)=rhoz;
GAM0(lambdastar,miu)=[rhomiu]*alpha/(1-alpha);
%===

% eq 9, capital utilization foc (Rk)
% -------------------------------------------
GAM0(Rk,Rk)=1;
GAM0(Rk,u)=-chi;

%===
GAM0(Rkstar,Rkstar)=1;
GAM0(Rkstar,ustar)=-chi;
%===

% eq 10, capital foc (phi)
% ------------------------------------------
GAM0(phi,phi)=1;
GAM0(phi,ephi)=-beta*exp(-gammastar-gammamiu)*(1-delta);
GAM0(phi,z)=rhoz;
GAM0(phi,miu)=[rhomiu]*[alpha/(1-alpha)+1];
GAM0(phi,elambda)=-(1-beta*exp(-gammastar-gammamiu)*(1-delta));
GAM0(phi,eRk)=-(1-beta*exp(-gammastar-gammamiu)*(1-delta));

%===
GAM0(phistar,phistar)=1;
GAM0(phistar,ephistar)=-beta*exp(-gammastar-gammamiu)*(1-delta);
GAM0(phistar,z)=rhoz;
GAM0(phistar,miu)=[rhomiu]*[alpha/(1-alpha)+1];
GAM0(phistar,elambdastar)=-(1-beta*exp(-gammastar-gammamiu)*(1-delta));
GAM0(phistar,eRkstar)=-(1-beta*exp(-gammastar-gammamiu)*(1-delta));
%===

% eq 11, investment foc (i)
% -------------------------------------------
expgmiu=exp(gammastar+gammamiu);
GAM0(i,lambda)=1/(S*expgmiu^2);
GAM0(i,phi)=-1/(S*expgmiu^2);
GAM0(i,upsilon)=-1/(S*expgmiu^2);
GAM0(i,miu)=((1-beta*rhomiu))*(alpha/(1-alpha)+1); 
GAM0(i,i)=(1+beta);
GAM0(i,z)=(1-beta*rhoz);
GAM0(i,ei)=-beta;
GAM1(i,i)=1;

%===
GAM0(istar,lambdastar)=1/(S*expgmiu^2);
GAM0(istar,phistar)=-1/(S*expgmiu^2);
GAM0(istar,upsilon)=-1/(S*expgmiu^2);
GAM0(istar,miu)=[(1-beta*rhomiu)]*[alpha/(1-alpha)+1]; 
GAM0(istar,istar)=(1+beta);
GAM0(istar,z)=(1-beta*rhoz);
GAM0(istar,eistar)=-beta;
GAM1(istar,istar)=1;
%===

% eq 12, capital utilization (k)
% -------------------------------------------
GAM0(k,k)=1;
GAM0(k,u)=-1;
GAM0(k,z)=1;
GAM1(k,kbar)=1;
GAM0(k,miu)=alpha/(1-alpha)+1;


%===
GAM0(kstar,kstar)=1;
GAM0(kstar,ustar)=-1;
GAM0(kstar,z)=1;
GAM1(kstar,kbarstar)=1;
GAM0(kstar,miu)=alpha/(1-alpha)+1;
%===

% eq 13, capital accumulation (kbar)
% ------------------------------------------
GAM0(kbar,kbar)=1;
GAM0(kbar,i)=-(1-(1-delta)*exp(-gammastar-gammamiu));
GAM0(kbar,upsilon)=-(1-(1-delta)*exp(-gammastar-gammamiu));
GAM1(kbar,kbar)=(1-delta)*exp(-gammastar-gammamiu);
GAM0(kbar,z)=(1-delta)*exp(-gammastar-gammamiu);
GAM0(kbar,miu)=((1-delta)*exp(-gammastar-gammamiu))*(alpha/(1-alpha)+1);

%===
GAM0(kbarstar,kbarstar)=1;
GAM0(kbarstar,istar)=-(1-(1-delta)*exp(-gammastar-gammamiu));
GAM0(kbarstar,upsilon)=-(1-(1-delta)*exp(-gammastar-gammamiu));
GAM1(kbarstar,kbarstar)=(1-delta)*exp(-gammastar-gammamiu);
GAM0(kbarstar,z)=(1-delta)*exp(-gammastar-gammamiu);
GAM0(kbarstar,miu)=((1-delta)*exp(-gammastar-gammamiu))*(alpha/(1-alpha)+1);
%===

% eq 14, goods market clearing (u)
% ------------------------------------------
GAM0(u,c)=css/yss;
GAM0(u,i)=iss/yss;
GAM0(u,y)=-1/gss;
GAM0(u,g)=1/gss;
GAM0(u,u)=kss*Rkss/yss;

%===
GAM0(ustar,cstar)=css/yss;
GAM0(ustar,istar)=iss/yss;
GAM0(ustar,ystar)=-1/gss;
GAM0(ustar,g)=1/gss;
GAM0(ustar,ustar)=kss*Rkss/yss;
%===

% eq 15, wage curve (wtilda)
% ------------------------------------------
AUX=(1+niu*(1+1/lambdaw))/(1-xiw*beta);
GAM0(wtilda,wtilda)=1;
GAM0(wtilda,x)=-1/AUX;
GAM0(wtilda,p)=iotaw*xiw*beta;
GAM0(wtilda,z)=-xiw*beta*(rhoz-iotaw);
GAM0(wtilda,miu)=[-xiw*beta*(rhomiu-iotaw)]*alpha/(1-alpha);
GAM0(wtilda,ewtilda)=-xiw*beta;
GAM0(wtilda,ep)=-xiw*beta;

% ===
AUXstar=(1+niu*(1+1/lambdaw));
GAM0(wtildastar,wtildastar)=1;
GAM0(wtildastar,xstar)=-1/AUXstar;
%===

% eq 16, MU labor (x)
% -----------------------------------------
GAM0(x,x)=1;
GAM0(x,psi)=-AUX*(1+beta)*xiw/(1-xiw);      
GAM0(x,b)=-1/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(x,L)=-niu;
GAM0(x,lambda)=1;
GAM0(x,w)=-niu*(1+1/lambdaw);

% ===
GAM0(xstar,xstar)=1;
GAM0(xstar,b)=-1/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(xstar,Lstar)=-niu;
GAM0(xstar,lambdastar)=1;
GAM0(xstar,wstar)=-niu*(1+1/lambdaw);
% ===

% eq 17, real wages indexation (w)
% ----------------------------------------
GAM0(w,w)=1;
GAM0(w,wtilda)=-(1-xiw);
GAM0(w,p)=xiw;
GAM0(w,z)=xiw;
GAM0(w,miu)=xiw*alpha/(1-alpha);
GAM1(w,w)=xiw;
GAM1(w,p)=iotaw*xiw;
GAM1(w,z)=iotaw*xiw;
GAM1(w,miu)=[iotaw*xiw]*alpha/(1-alpha);

% ===
GAM0(wstar,wstar)=1;
GAM0(wstar,wtildastar)=-1;
% ===

% eq 18, MP rule (R)
% ----------------------------------------
GAM0(R,R)=1;
GAM1(R,R)=rhoR;
GAM0(R,p)=-(1-rhoR)*fp;
GAM0(R,gdp)=-(1-rhoR)*fy-fdy;
GAM0(R,gdpstar)=(1-rhoR)*fy+fdy;
GAM1(R,gdp)=-fdy;
GAM1(R,gdpstar)=fdy;
GAM0(R,mp)=-1;

% eq 19 - 24, exogenous shocks 
% ------------------------------------------
GAM0(z,z)=1; GAM1(z,z)=rhoz; PSI(z,zs)=1;
GAM0(g,g)=1; GAM1(g,g)=rhog; PSI(g,gs)=1;

GAM0(psi,psi)=1; GAM1(psi,psi)=rhopsi; GAM0(psi,ARMApsi)=-1; GAM1(psi,ARMApsi)=-rhoARMApsi;
GAM0(ARMApsi,ARMApsi)=1; PSI(ARMApsi,psis)=1;

GAM0(miu,miu)=1; GAM1(miu,miu)=rhomiu; PSI(miu,mius)=1;

GAM0(lambdap,lambdap)=1; GAM1(lambdap,lambdap)=rholambdap; GAM0(lambdap,ARMAlambdap)=-1; GAM1(lambdap,ARMAlambdap)=-rhoARMAlambdap;
GAM0(ARMAlambdap,ARMAlambdap)=1; PSI(ARMAlambdap,lambdaps)=1;

GAM0(b,b)=1; GAM1(b,b)=rhob; PSI(b,bs)=1;
GAM0(mp,mp)=1; GAM1(mp,mp)=rhomp; PSI(mp,Rs)=1;
GAM0(upsilon,upsilon)=1; GAM1(upsilon,upsilon)=rhoupsilon; PSI(upsilon,upsilons)=1;

% eq 25 - 32, expectational terms
% ------------------------------------------
GAM0(ep,p)=1; GAM1(ep,ep)=1; PPI(ep,pex)=1;
GAM0(ec,c)=1; GAM1(ec,ec)=1; PPI(ec,cex)=1;
GAM0(elambda,lambda)=1; GAM1(elambda,elambda)=1; PPI(elambda,lambdaex)=1;
GAM0(ephi,phi)=1; GAM1(ephi,ephi)=1; PPI(ephi,phiex)=1;
GAM0(eRk,Rk)=1; GAM1(eRk,eRk)=1; PPI(eRk,Rkex)=1;
GAM0(ei,i)=1; GAM1(ei,ei)=1; PPI(ei,iex)=1;
GAM0(ewtilda,wtilda)=1; GAM1(ewtilda,ewtilda)=1; PPI(ewtilda,wtildaex)=1;

%===
GAM0(ecstar,cstar)=1; GAM1(ecstar,ecstar)=1; PPI(ecstar,cexstar)=1;
GAM0(elambdastar,lambdastar)=1; GAM1(elambdastar,elambdastar)=1; PPI(elambdastar,lambdaexstar)=1;
GAM0(ephistar,phistar)=1; GAM1(ephistar,ephistar)=1; PPI(ephistar,phiexstar)=1;
GAM0(eRkstar,Rkstar)=1; GAM1(eRkstar,eRkstar)=1; PPI(eRkstar,Rkexstar)=1;
GAM0(eistar,istar)=1; GAM1(eistar,eistar)=1; PPI(eistar,iexstar)=1;
%===

% Observables
% -------------
GAM0(yobs,yobs)=1; 
GAM0(yobs,gdp)= -1;
GAM0(yobs,z)= -1; 
GAM0(yobs,miu)= -alpha/(1-alpha); 
GAM1(yobs,gdp)= -1; 

GAM0(cobs,cobs)=1; 
GAM0(cobs,c)= -1;
GAM0(cobs,z)= -1; 
GAM0(cobs,miu)= -alpha/(1-alpha); 
GAM1(cobs,c)= -1; 

GAM0(iobs,iobs)=1; 
GAM0(iobs,i)= -1;
GAM0(iobs,z)= -1; 
GAM0(iobs,miu)= -alpha/(1-alpha); 
GAM1(iobs,i)= -1; 

GAM0(wobs,wobs)=1; 
GAM0(wobs,w)= -1;
GAM0(wobs,z)= -1; 
GAM0(wobs,miu)= -alpha/(1-alpha); 
GAM1(wobs,w)= -1; 

GAM0(dagtrend,dagtrend)= -1; 
GAM0(dagtrend,z)=1; 
GAM0(dagtrend,miu)=alpha/(1-alpha); 

GAM0(ygaplev,ygaplev)= -1; 
GAM0(ygaplev,gdp)=1;
GAM0(ygaplev,gdpstar)= -1; 

% Standard Deviations 
if flag_sv == 0 
    SDXone=diag(param(ncof+1:end)); 
else 
    SDXone=[]; 
end 
Cone=zeros(NY,1) ;

[Gone,~,RRone,~,~,~,gev,euone]=gensys(GAM0,GAM1,Cone,PSI,PPI);
if ~isequal(euone,[1;1])
    Gtwo=[]; Ctwo=[]; RRtwo=[]; eutwo=[0;0]; SDXtwo=[]; zmat=[]; 
    gev=[]; Cone=[]; 
    return; 
end 


% ZMAT matrix maps observables to states 
% Size depends on whether relative price observable or not 
if flag_relpiobs==1
    zmat=zeros(8,size(Gone,2));
else
    zmat=zeros(7,size(Gone,2));
end
zmat(1,yobs)=1; 
zmat(2,cobs)=1; 
zmat(3,iobs)=1; 
zmat(4,L)=1; 
zmat(5,wobs)=1; 
zmat(6,p)=1;
zmat(7,R)=1;

% Vector of constants in observation equation 
%==========================
Cone=zeros(NY,1); 
Cone(yobs)=gammastar100;
Cone(cobs)=gammastar100;
Cone(iobs)=gammastar100;
Cone(L)=Lss;
Cone(wobs)=gammastar100;
Cone(p)=pss100;
Cone(R)=pss100+rss100;

if flag_relpiobs==1 
    zmat(8,miu)=1; 
    Cone(miu)=gammamiu100; 
end 

% =======================
% Solution, second sample 
% =======================
[Gtwo,Ctwo,RRtwo,eutwo,SDXtwo]=modptjREDsub(partwo);
if ~isequal(eutwo,[1;1])
	Gtwo=[];zmat=[];return; 
end 

% Verified that if make coefficient above equal across subsamples solutions
% are the same Jan 30 2009
% comparemat(Gone,Gtwo); 
% comparemat(RRone,RRtwo); 
% comparemat(zmat*Cone,zmat2*Ctwo)

if ~isequal(size(SDXone),size(SDXtwo)); 
	error('Size of Covariance Matrices do not match') 
end 

ncof=ncof+1; 