@#include "model_declarations.mod"

//----------------------------------------------------------------
// 4. Initial values for parameters
//----------------------------------------------------------------
bbeta        = 0.99;
kkappa       = 0.50;
inv_ssigma   = 1.00;
rrho_i       = 0.85;
pphi_pi      = 1.50;
pphi_ytilde  = 0.5/4;

 
   
    

ppsi_yan     = 1;

rrho_v       = 0.0;
rrho_a       = 0.9;
rrho_u       = 0.9;


set_param_value('ssigma_v',0.1/100);
set_param_value('ssigma_a',0.05/100);
set_param_value('ssigma_u',0.05/100);

 

rrho         = -log(bbeta);


@#include "model_equations.mod"    

resid;
check;
steady;


shocks;
    var vareps_v;
     stderr 1;
    var vareps_a;
     stderr 1;
    var vareps_u;
     stderr 1;
end;



set_dynare_seed(unifrnd(0,1e8));
stoch_simul(order=1,periods=15,drop=0,irf=0,noprint,nograph) ytildeobs ppiobs Robs;




