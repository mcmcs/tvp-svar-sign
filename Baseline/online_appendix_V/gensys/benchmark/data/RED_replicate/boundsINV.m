function theta = boundsINV (param)

alpha=bound01INV(param(1));     
iotap=bound01INV(param(2));     
iotaw=bound01INV(param(3));     
gammastar100=param(4);  
gammamiu100 =param(5); 
h=bound01INV(param(6));  
lambdapss=param(7); 
lambdawss=param(8); 
Lss=param(9);       
pss100=param(10);   

Fbeta=bound0INV(param(11));
niu=bound0INV(param(12));
xip=bound01INV(param(13));
xiw=bound01INV(param(14));
chi=bound0INV(param(15));
S=bound0INV(param(16));
fp=param(17); 
fy=param(18);

rhoR=modtomin_ab(param(19),0,.99);
rhoz=modtomin_ab(param(20),0,.99);
rhog=modtomin_ab(param(21),0,.99);
rhomiu=modtomin_ab(param(22),0,.99);
rholambdap=modtomin_ab(param(23),0,.99);
rholambdaw=modtomin_ab(param(24),0,.99);
rhob=modtomin_ab(param(25),0,.99);
rhoARMAlambdap=modtomin_ab(param(26),0,.99);
rhoARMAlambdaw=modtomin_ab(param(27),0,.99);
fdy=param(28); 
rhoMEI=modtomin_ab(param(29),0,.99); 

gammastars2=param(30); 
gammamius2 =param(31); 

sdR=bound0INV(param(32));
sdz=bound0INV(param(33)); 
sdg=bound0INV(param(34)); 
sdmiu=bound0INV(param(35));
sdlambdap=bound0INV(param(36)); 
sdlambdaw=bound0INV(param(37)); 
sdb=bound0INV(param(38));
sdMEI=bound0INV(param(39));

theta = [alpha iotap iotaw gammastar100 gammamiu100 h lambdapss lambdawss Lss pss100...
    Fbeta niu xip xiw chi S fp fy  rhoR rhoz...
    rhog rhomiu rholambdap rholambdaw rhob rhoARMAlambdap rhoARMAlambdaw fdy rhoMEI gammastars2...
    gammamius2 sdR sdz sdg sdmiu sdlambdap sdlambdaw sdb sdMEI]';

function rho = bound01INV(param);
rho = log(param/(1-param));

function sigma = bound0INV(param);
sigma = log(param);