@#include "model_declarations.mod"

//----------------------------------------------------------------
// 4. Initial values for parameters
//----------------------------------------------------------------
bbeta        = 0.99;
kkappa       = 0.5;
inv_ssigma   = 1.00;
rrho_i       = 0.85;
pphi_pi      = 1.5;
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







//varobs ytildeobs Robs ppiobs; 
resid;
check;



stoch_simul(order=1,irf=20) ytilde ppi R;

switch kkappa

case 0.1
oo_flat = oo_;

save('flat_PC.mat','oo_flat')
case 0.5

oo_steep = oo_;

save('steep_PC.mat','oo_steep')
 end