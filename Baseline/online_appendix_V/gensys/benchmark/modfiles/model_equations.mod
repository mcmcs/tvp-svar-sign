model; 


//1. Euler equation
ytilde = -inv_ssigma*(R-ppi(+1)-rn) + ytilde(+1);

//2. PC
ppi = bbeta*ppi(+1) + kkappa*ytilde + u;

//3. MP equation
R = rrho_i*R(-1) + (1-rrho_i)*(rrho + pphi_pi*ppi +pphi_ytilde*ytilde) + v;

//4. Natural rate of interest
rn = rrho + (1/inv_ssigma)*ppsi_yan*(rrho_a-1)*a;

//5. Monetary policy shock
v=rrho_v*v(-1)  + ssigma_v*vareps_v;

//6. Technology shock
a=rrho_a*a(-1)  + ssigma_a*vareps_a;

//7. Cost push shock
u=rrho_u*u(-1)  + ssigma_u*vareps_u;



end;


initval;
    ytilde    = 0;
    ppi       = 0;
    R         = rrho;
    rn        = rrho;
    v         = 0;
    a         = 0;
    u         = 0;
end;


