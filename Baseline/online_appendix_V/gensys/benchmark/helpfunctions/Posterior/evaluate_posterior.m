function f = evaluate_posterior(ttheta,Y,pmean, pstdd, pshape)


% check that the parameters satisfy are in the domain if not, set f=-1e8 and
% and draw a new set of parameters


ret = get_parameter_bounds(ttheta);


if ret<1
    f=-1e8;
    
else
    
    [GAM0, GAM1, PPSI,PPI,C] = get_model(ttheta);           % Linearized model evaluated at ttheta0
    [T1,~,T0,~,~,~,~,eu,~]   = gensys(GAM0,GAM1,C,PPSI,PPI); % Policy functions evaluated at ttheta0
    
    
    
    % record determinacy | indeterminacy | non-existence
    if (eu(1)==1 && eu(2)==1)
        
        
        loglh   = evaluate_like(ttheta,Y,T1,T0);
        
        lnprior = evaluate_prior(ttheta, pmean, pstdd, pshape);
        
        
        
    else
        
        
        % there is no equilibrium or equilibrium is not unique
        loglh   = -1e8;
        lnprior = -1e8;
    end
    
    
    % checking prior density
    f = loglh + lnprior;
    
end
