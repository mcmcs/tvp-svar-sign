function theta = bounds(param)
% Take vector from constrained minimization and tranform it into vector for
% model solution. Only estimated parameters. 
% Inverse function is boundsINV.m 
alpha=bound01(param(1));        
iotap=bound01(param(2));        
iotaw=bound01(param(3));        
gammastar100=param(4);          
gammamiu100 =param(5);          
h=bound01(param(6));  
lambdapss=param(7); 
lambdawss=param(8); 
Lss=param(9);       
pss100=param(10);

Fbeta=bound0(param(11));
niu=bound0(param(12));
xip=bound01(param(13));
xiw=bound01(param(14));
chi=bound0(param(15));
S=bound0(param(16));
fp=param(17); 
fy=param(18);       

rhoR=mintomod_ab(param(19),0,.99);
rhoz=mintomod_ab(param(20),0,.99);
rhog=mintomod_ab(param(21),0,.99);
rhomiu=mintomod_ab(param(22),0,.99);
rholambdap=mintomod_ab(param(23),0,.99);
rholambdaw=mintomod_ab(param(24),0,.99);
rhob=mintomod_ab(param(25),0,.99);
rhoARMAlambdap=mintomod_ab(param(26),0,.99);
rhoARMAlambdaw=mintomod_ab(param(27),0,.99);
fdy=param(28);      
rhoMEI=mintomod_ab(param(29),0,.99); 

gammastars2=param(30); 
gammamius2 =param(31); 

sdR=bound0(param(32));
sdz=bound0(param(33)); 
sdg=bound0(param(34)); 
sdmiu=bound0(param(35));
sdlambdap=bound0(param(36)); 
sdlambdaw=bound0(param(37)); 
sdb=bound0(param(38));
sdMEI=bound0(param(39)); 

theta = [alpha iotap iotaw gammastar100 gammamiu100 h lambdapss lambdawss Lss pss100...
    Fbeta niu xip xiw chi S fp fy  rhoR rhoz...
    rhog rhomiu rholambdap rholambdaw rhob rhoARMAlambdap rhoARMAlambdaw fdy rhoMEI gammastars2...
    gammamius2 sdR sdz sdg sdmiu sdlambdap sdlambdaw sdb sdMEI]';

function rho = bound01(param);
rho = 1-1/(1+exp(param));


function sigma = bound0(param);
sigma = exp(param);