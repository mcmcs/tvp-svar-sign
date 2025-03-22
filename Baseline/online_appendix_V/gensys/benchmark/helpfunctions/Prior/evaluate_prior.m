function lnprior = evaluate_prior(para, pmean, pstdd, pshape)


nprio    = size(pshape,1);
prioinfo = horzcat(zeros(nprio,2),pshape);



lnprior        = 0;

i = 1;
while i <= nprio
    
    switch prioinfo(i,3)
        
        case 1 % BETA Prior
            
            
            a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
            b = a*(1/pmean(i) - 1);
            
            lnprior = lnprior + lpdfbeta(para(i),a,b);

            
        case 2 % GAMMA PRIOR
            
            b = pstdd(i)^2/pmean(i);
            a = pmean(i)/b;
            
            
            lnprior = lnprior + lpdfgam(para(i),a,b);

            
        case 3 % GAUSSIAN PRIOR
            
            a = pmean(i);
            b = pstdd(i);
            
            lnprior = lnprior + lpdfnor(para(i),a,b);

            
        case 4 % INVGAMMA PRIOR
            
            a = 2 + (pmean(i)/pstdd(i))^2;
            b = (a-1)*pmean(i);
            
            lnprior = lnprior + lpdfig(para(i),a,b);
            
                        
            
    end;
    
    
    i = i+1;
    
    
end

end


