function [lnprior,lnprior_matlab,my_densities,matlab_densities] = get_prior_density_check(para, pmean, pstdd, pshape)


nprio    = size(pshape,1);
prioinfo = horzcat(zeros(nprio,2),pshape);

my_densities     = zeros(nprio,1);
matlab_densities = zeros(nprio,1);

lnprior        = 0;
lnprior_matlab = 0;


i = 1;
while i <= nprio
    
    
    
    
    switch prioinfo(i,3)
        
        case 1 % BETA Prior
            
            
            a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
            b = a*(1/pmean(i) - 1);
            
            lnprior = lnprior + lpdfbeta(para(i),a,b);
            lnprior_matlab = lnprior_matlab + log(betapdf(para(i),a,b));
            
            
            my_densities(i,1)      = lpdfbeta(para(i),a,b);
            matlab_densities(i,1)  = log(betapdf(para(i),a,b));
            
        case 2 % GAMMA PRIOR
            
            b = pstdd(i)^2/pmean(i);
            a = pmean(i)/b;
            
            
            lnprior = lnprior + lpdfgam(para(i),a,b);
            lnprior_matlab = lnprior_matlab + lpdfgam(para(i),a,b);
            
            my_densities(i,1)      = lpdfgam(para(i),a,b);
            matlab_densities(i,1)  = log(gampdf(para(i),a,b));
            
        case 3 % GAUSSIAN PRIOR
            
            a = pmean(i);
            b = pstdd(i);
            
            lnprior = lnprior + lpdfnor(para(i),a,b);
            lnprior_matlab = lnprior_matlab + log(normpdf(para(i),a,b));
            
            my_densities(i,1)      = lpdfnor(para(i),a,b);
            matlab_densities(i,1)  = log(normpdf(para(i),a,b));
            
        case 4 % INVGAMMA PRIOR
            
            a = 2 + (pmean(i)/pstdd(i))^2;
            b = (a-1)*pmean(i);
            
            lnprior = lnprior + lpdfig(para(i),a,b);
            
            jac                    = log((1/para(i))^2);
            
            lnprior_matlab         = lnprior_matlab + (log(gampdf(1/para(i),a,1/b)) +jac);
            
            my_densities(i,1)      = lpdfig(para(i),a,b);
            
            matlab_densities(i,1)  = (log(gampdf(1/para(i),a,1/b)) + jac);
            
            
            
    end;
    
    prioinfo(i,1) = a;
    prioinfo(i,2) = b;
    
    
    i = i+1;
    
    
end

end


