function ret = get_parameter_bounds(ttheta)

domain = [1e-5 20;                % 1 psippi
    1e-5 20;                      % 2 ppsiy
    1e-5 20;                      % 3 ppsigy
    -2 2;                         % 4 rrhoR1
    -2 2;                         % 5 rrhoR2
    1e-5 20;                      % 6 kkappa
    1e-5 0.99999;                 % 7 b
    1e-5   1;                     % 8 eeta
    1e-5   1;                     % 9 eetaw
    1e-5 20;                      % 10 ttau
    1e-5 0.99999;                 % 11 nnu
    1e-5 0.99999;                 % 12 nnuw
    1e-5 0.99999;                 % 13 rrrhoG
    1e-5 0.99999;                 % 14 rrhoZ
    1e-5 0.99999;                 % 15 rrhod
    1e-5 0.99999;                 % 16 rrhodL
    1e-5 0.99999;                 % 17 rrhomI
    1e-5 100;                     % 18 ssigmad
    1e-7 100;                     % 19 ssigmadL
    1e-7 100;                     % 20 ssigmaZ
    1e-7 100;                     % 21 ssigmaMI
    1e-7 100;                     % 22 ssigmaR
    1e-7 100;                     % 23 ssigamG
    1e-7 100;                     % 24 ssigmaobsppi
    1e-7 100;                     % 25 ssigmaobsw
    ];


para_bd1 = ttheta>=domain(:,1);
para_bd2 = ttheta<=domain(:,2);


% if any parameter restriction is violated then ret <1
ret = (sum(para_bd1) + sum(para_bd2))/(size(ttheta,1)*2);

if  abs(ttheta(4) + ttheta(5))  > 1
    ret = 0;
end

end
      
