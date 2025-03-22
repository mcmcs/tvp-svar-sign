function out = plot_prior_draws(prior_draws,pmean, pstdd, pshape)

close all

nprio    = size(pshape,1);
prioinfo = horzcat(zeros(nprio,2),pshape);

%% Monetary Policy Equation Parameters
break1 = 5;
hFig_MP = figure('name','Monetary Policy Equation Parameters','NumberTitle','off');
set(hFig_MP, 'Position', [0 20 650 325])
for i=1:break1
    
    subplot(2,3,i)
    H = histogram(prior_draws(i,:),'normalization','probability');
    nbinedges = length(H.BinEdges);
    
    densities = zeros(size(prior_draws,1),nbinedges);
    
    
    
    for j=1:nbinedges
        
        
        
        switch prioinfo(i,3)
            
            case 1 % BETA Prior
                
                a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
                b = a*(1/pmean(i) - 1);
                
                densities(i,j)      = exp(lpdfbeta(H.BinEdges(j),a,b));
                
                
            case 2 % GAMMA PRIOR
                
                b = pstdd(i)^2/pmean(i);
                a = pmean(i)/b;
                
                densities(i,j)      = exp(lpdfgam(H.BinEdges(j),a,b));
                
                
            case 3 % GAUSSIAN PRIOR
                
                a = pmean(i);
                b = pstdd(i);
                
                densities(i,j)      = exp(lpdfnor(H.BinEdges(j),a,b));
                
            case 4 % INVGAMMA PRIOR
                
                a = 2 + (pmean(i)/pstdd(i))^2;
                b = (a-1)*pmean(i);
                
                % draw from Gamma(gam_a,gam_b)
                gam_a = a;
                gam_b = 1/b;
                
                
                %jac                    = log((1/prior_draws(i,j))^2);
                
                densities(i,j)      = exp(lpdfig(H.BinEdges(j),a,b));
                
                % matlab_densities(i,j)  = (log(gampdf(1/prior_draws(i,j),a,1/b)) + jac);
                
                
                
        end
        
        
        
    end
    
    
    hold on
    plot(H.BinEdges,(densities(i,:)/max(densities(i,:)))*max(H.Values))
    
    switch i
        case 1
            
            title('$\psi_{\pi}$','interpreter','latex')
            
        case 2
            
            title('$\psi_{y}$','interpreter','latex')
            
        case 3
            
            title('$\psi_{gy}$','interpreter','latex')
        case 4
            
            title('$\rho_{1}$','interpreter','latex')
            
        case 5
            
            title('$\rho_{2}$','interpreter','latex')
            
        otherwise
    end
    
    
end





%% Structural parameters
break2 = break1+8;
hFig_structural_para = figure('name','Structural parameters','NumberTitle','off');
set(hFig_structural_para, 'Position', [0 20 550 725])
for i=(break1+1):break2
    
    subplot(4,2,i-break1)
    H = histogram(prior_draws(i,:),'normalization','probability');
    nbinedges = length(H.BinEdges);
    
    densities = zeros(size(prior_draws,1),nbinedges);
    
    
    
    for j=1:nbinedges
        
        
        
        switch prioinfo(i,3)
            
            case 1 % BETA Prior
                
                a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
                b = a*(1/pmean(i) - 1);
                
                densities(i,j)      = exp(lpdfbeta(H.BinEdges(j),a,b));
                
                
            case 2 % GAMMA PRIOR
                
                b = pstdd(i)^2/pmean(i);
                a = pmean(i)/b;
                
                densities(i,j)      = exp(lpdfgam(H.BinEdges(j),a,b));
                
                
            case 3 % GAUSSIAN PRIOR
                
                a = pmean(i);
                b = pstdd(i);
                
                densities(i,j)      = exp(lpdfnor(H.BinEdges(j),a,b));
                
            case 4 % INVGAMMA PRIOR
                
                a = 2 + (pmean(i)/pstdd(i))^2;
                b = (a-1)*pmean(i);
                
                % draw from Gamma(gam_a,gam_b)
                gam_a = a;
                gam_b = 1/b;
                
                
                %jac                    = log((1/prior_draws(i,j))^2);
                
                densities(i,j)      = exp(lpdfig(H.BinEdges(j),a,b));
                
                % matlab_densities(i,j)  = (log(gampdf(1/prior_draws(i,j),a,1/b)) + jac);
                
                
                
        end
        
        
        
    end
    
    
    hold on
    plot(H.BinEdges,(densities(i,:)/max(densities(i,:)))*max(H.Values))
    
    switch i
        case 6
            
            title('$\kappa$','interpreter','latex')
            
        case 7
            
            title('$b$','interpreter','latex')
            
        case 8
            
            title('$\eta_{p}$','interpreter','latex')
        case 9
            
            title('$\eta_{w}$','interpreter','latex')
            
        case 10
            
            title('$\tau$','interpreter','latex')
        case 11
            
            title('$\beta$','interpreter','latex')
            
        case 12
            
            title('$\nu_{p}$','interpreter','latex')
            
        case 13
            
            title('$\nu_{w}$','interpreter','latex')
            
        otherwise
    end
    
    
end







%% Exogenous processes: persistance
break3 = break2+6;
hFig_exo_persistance = figure('name','Exogenous processes: persistance','NumberTitle','off');
set(hFig_exo_persistance, 'Position', [0 20 650 325])
for i=(break2+1):break3
    
    subplot(2,3,i-break2)
    H = histogram(prior_draws(i,:),'normalization','probability');
    nbinedges = length(H.BinEdges);
    
    densities = zeros(size(prior_draws,1),nbinedges);
    
    
    
    for j=1:nbinedges
        
        
        
        switch prioinfo(i,3)
            
            case 1 % BETA Prior
                
                a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
                b = a*(1/pmean(i) - 1);
                
                densities(i,j)      = exp(lpdfbeta(H.BinEdges(j),a,b));
                
                
            case 2 % GAMMA PRIOR
                
                b = pstdd(i)^2/pmean(i);
                a = pmean(i)/b;
                
                densities(i,j)      = exp(lpdfgam(H.BinEdges(j),a,b));
                
                
            case 3 % GAUSSIAN PRIOR
                
                a = pmean(i);
                b = pstdd(i);
                
                densities(i,j)      = exp(lpdfnor(H.BinEdges(j),a,b));
                
            case 4 % INVGAMMA PRIOR
                
                a = 2 + (pmean(i)/pstdd(i))^2;
                b = (a-1)*pmean(i);
                
                % draw from Gamma(gam_a,gam_b)
                gam_a = a;
                gam_b = 1/b;
                
                
                %jac                    = log((1/prior_draws(i,j))^2);
                
                densities(i,j)      = exp(lpdfig(H.BinEdges(j),a,b));
                
                % matlab_densities(i,j)  = (log(gampdf(1/prior_draws(i,j),a,1/b)) + jac);
                
                
                
        end
        
        
        
    end
    
    
    hold on
    plot(H.BinEdges,(densities(i,:)/max(densities(i,:)))*max(H.Values))
    
    switch i
        case 14
            
            title('$\rho_{d}$','interpreter','latex')
            
        case 15
            
            title('$\rho_{d_{L}}$','interpreter','latex')
            
        case 16
            
            title('$\rho_{A}$','interpreter','latex')
        case 17
            
            title('$\rho_{m_{I}}$','interpreter','latex')
            
        case 18
            
            title('$\rho_{\mu}$','interpreter','latex')
        case 19
            
            title('$\rho_{G}$','interpreter','latex')
            
            
        otherwise
    end
    
    
end





%% Exogenous processes: standard deviation
break4 = break3 + 7;
hFig_exo_std = figure('name','Exogenous processes: Standard Deviation','NumberTitle','off');
set(hFig_exo_std , 'Position', [0 20 550 725])
for i=(break3+1):break4
    
    subplot(4,2,i-break3)
    H = histogram(prior_draws(i,:),'normalization','probability');
    nbinedges = length(H.BinEdges);
    
    densities = zeros(size(prior_draws,1),nbinedges);
    
    
    
    for j=1:nbinedges
        
        
        
        switch prioinfo(i,3)
            
            case 1 % BETA Prior
                
                a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
                b = a*(1/pmean(i) - 1);
                
                densities(i,j)      = exp(lpdfbeta(H.BinEdges(j),a,b));
                
                
            case 2 % GAMMA PRIOR
                
                b = pstdd(i)^2/pmean(i);
                a = pmean(i)/b;
                
                densities(i,j)      = exp(lpdfgam(H.BinEdges(j),a,b));
                
                
            case 3 % GAUSSIAN PRIOR
                
                a = pmean(i);
                b = pstdd(i);
                
                densities(i,j)      = exp(lpdfnor(H.BinEdges(j),a,b));
                
            case 4 % INVGAMMA PRIOR
                
                a = 2 + (pmean(i)/pstdd(i))^2;
                b = (a-1)*pmean(i);
                
                % draw from Gamma(gam_a,gam_b)
                %gam_a = a;
                %gam_b = 1/b;
                %jac                    = log((1/prior_draws(i,j))^2);
                
                densities(i,j)      = exp(lpdfig(H.BinEdges(j),a,b));
                
                % matlab_densities(i,j)  = (log(gampdf(1/prior_draws(i,j),a,1/b)) + jac);
                
                
                
        end
        
        
        
    end
    
    
    hold on
    plot(H.BinEdges,(densities(i,:)/max(densities(i,:)))*max(H.Values))
    
    switch i
        case 20
            
            title('$\sigma_{d}$','interpreter','latex')
            
        case 21
            
            title('$\sigma_{d_{L}}$','interpreter','latex')
            
        case 22
            
            title('$\sigma_{A}$','interpreter','latex')
        case 23
            
            title('$\sigma_{m_{I}}$','interpreter','latex')
            
        case 24
            
            title('$\sigma_{\mu}$','interpreter','latex')

        
        case 25
            
            title('$\sigma_{R}$','interpreter','latex')   
            
        case 26
            
            title('$\sigma_{G}$','interpreter','latex')
            
            
        otherwise
    end
    
    
end






%% Measurement equation shocks
break5 = break4 + 2;
hFig_exo_std = figure('name','Measurement equation: Standard Deviation','NumberTitle','off');
set(hFig_exo_std , 'Position', [0 20 450 225])
for i=(break4+1):break5
    
    subplot(1,2,i-break4)
    H = histogram(prior_draws(i,:),'normalization','probability');
    nbinedges = length(H.BinEdges);
    
    densities = zeros(size(prior_draws,1),nbinedges);
    
    
    
    for j=1:nbinedges
        
        
        
        switch prioinfo(i,3)
            
            case 1 % BETA Prior
                
                a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
                b = a*(1/pmean(i) - 1);
                
                densities(i,j)      = exp(lpdfbeta(H.BinEdges(j),a,b));
                
                
            case 2 % GAMMA PRIOR
                
                b = pstdd(i)^2/pmean(i);
                a = pmean(i)/b;
                
                densities(i,j)      = exp(lpdfgam(H.BinEdges(j),a,b));
                
                
            case 3 % GAUSSIAN PRIOR
                
                a = pmean(i);
                b = pstdd(i);
                
                densities(i,j)      = exp(lpdfnor(H.BinEdges(j),a,b));
                
            case 4 % INVGAMMA PRIOR
                
                a = 2 + (pmean(i)/pstdd(i))^2;
                b = (a-1)*pmean(i);
                
                % draw from Gamma(gam_a,gam_b)
                gam_a = a;
                gam_b = 1/b;
                
                
                %jac                    = log((1/prior_draws(i,j))^2);
                
                densities(i,j)      = exp(lpdfig(H.BinEdges(j),a,b));
                
                % matlab_densities(i,j)  = (log(gampdf(1/prior_draws(i,j),a,1/b)) + jac);
                
                
                
        end
        
        
        
    end
    
    
    hold on
    plot(H.BinEdges,(densities(i,:)/max(densities(i,:)))*max(H.Values))
    
    switch i
        case 27
            
            title('$\sigma_{obs,\pi}$','interpreter','latex')
            
        case 28
            
            title('$\sigma_{obs,w}$','interpreter','latex')

            
        otherwise
    end
    
    
end

out =1;

end

