function [prior_draws,my_densities,matlab_densities] = get_prior_draws(ndraws,pmean, pstdd, pshape)


nprio    = size(pshape,1);
prioinfo = horzcat(zeros(nprio,2),pshape);

prior_draws      = zeros(nprio,ndraws);

my_densities     = zeros(nprio,ndraws);
matlab_densities = zeros(nprio,ndraws);





for j=1:ndraws
    
    i = 1;
    while i <= nprio
        
        
        
        
        switch prioinfo(i,3)
            
            case 1 % BETA Prior
                
                
                a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
                b = a*(1/pmean(i) - 1);
                
                prior_draws(i,j) = betarnd(a,b);
                
                my_densities(i,j)      = lpdfbeta(prior_draws(i,j),a,b);
                matlab_densities(i,j)  = log(betapdf(prior_draws(i,j),a,b));
                
            case 2 % GAMMA PRIOR
                
                b = pstdd(i)^2/pmean(i);
                a = pmean(i)/b;
                
                prior_draws(i,j) = gamrnd(a,b);
                
                my_densities(i,j)      = lpdfgam(prior_draws(i,j),a,b);
                matlab_densities(i,j)  = log(gampdf(prior_draws(i,j),a,b));
                
            case 3 % GAUSSIAN PRIOR
                
                a = pmean(i);
                b = pstdd(i);
                
                prior_draws(i,j) = normrnd(a,b);
                
                my_densities(i,j)      = lpdfnor(prior_draws(i,j),a,b);
                matlab_densities(i,j)  = log(normpdf(prior_draws(i,j),a,b));
                
            case 4 % INVGAMMA PRIOR
                
                a = 2 + (pmean(i)/pstdd(i))^2;
                b = (a-1)*pmean(i);
                
                % draw from Gamma(gam_a,gam_b)
                gam_a = a;
                gam_b = 1/b;
                
                prior_draws(i,j) = 1/gamrnd(gam_a,gam_b);
                
                
                jac                    = log((1/prior_draws(i,j))^2);
                
                my_densities(i,j)      = lpdfig(prior_draws(i,j),a,b);
                
                matlab_densities(i,j)  = (log(gampdf(1/prior_draws(i,j),a,1/b)) + jac);
                
                
                
        end;
        
        
        i = i+1;
        
        
    end
    
end


end


