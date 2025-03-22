function [ttheta0,H0] = get_initial_guess(initial_guess)



if initial_guess==1
    
    try
        
        %load '~/Dropbox/trendinflation_2017/code/benchmark/matfiles/optimizedttheta.mat';
        load('~/Dropbox/trendinflation_2017/code/benchmark/estimation_dynare/my_driver_results.mat');
        
        ttheta0    = ...
        [oo_.posterior_mode.parameters.ppsipi          
        oo_.posterior_mode.parameters.ppsiy            
     oo_.posterior_mode.parameters.ppsigy           
     oo_.posterior_mode.parameters.rrhoR1           
     oo_.posterior_mode.parameters.rrhoR2           
     oo_.posterior_mode.parameters.kkappa           
     oo_.posterior_mode.parameters.b                
     oo_.posterior_mode.parameters.transeeta     
     oo_.posterior_mode.parameters.transeetaw       
     oo_.posterior_mode.parameters.ttau                
     oo_.posterior_mode.parameters.nnu              
     oo_.posterior_mode.parameters.nnuw             
     oo_.posterior_mode.parameters.rrhod            
     oo_.posterior_mode.parameters.rrhodL           
     oo_.posterior_mode.parameters.rrhoA           
     oo_.posterior_mode.parameters.rrhomI                
     oo_.posterior_mode.parameters.rrhoG            
     oo_.posterior_mode.shocks_std.vareepsilond    
     oo_.posterior_mode.shocks_std.vareepsilondL  
     oo_.posterior_mode.shocks_std.vareepsilonA     
     oo_.posterior_mode.shocks_std.vareepsilonmI      
     oo_.posterior_mode.shocks_std.vareepsilonR 
     oo_.posterior_mode.shocks_std.vareepsilonG 
     oo_.posterior_mode.shocks_std.obserrorppi 
     oo_.posterior_mode.shocks_std.obserrorw    ];
 

        
    catch
        
        display('There does not exist and optimized guess. Please set initial_guess = 0')
        
    end
    
else



    ttheta0=[
2.279749361681230; 
0.013752272023865;
0.524235424966235;
1.290185481537505;
-0.434799703883344;
3.982292138265899;
0.820917431257614;
0.261972418750297;
0.204500139701313;
1.134147303218401;
0.816418765319189;
0.575029895195554;
0.480129369809981;
0.967915284467965;
0.095369575530322;
0.710606308501769;
0.965756204932250;
0.020169461976973;
0.020845285511036;
0.005295949101297;
0.023811203386023;
0.0003;
0.001019542013110;
0.016467102351049;
0.001443932117735;
0.006826136282733];    % para 25 std shock measurement equation for real wage growth

end


    c0=1e-3;
    npara = size(ttheta0,1);
    H0=c0*eye(npara);





end


