function prior = get_prior_specification()

prior = [
1.7000  0.3000  2.0000   % para 1  ppsi_ppi response to inflation 
0.1300  0.1000  2.0000   % para 2  ppsi_gy  response to output growth
0.1300  0.1000  2.0000   % para 3  ppsi_y
1.0000  0.2000  3.0000   % para 4  rrho_R1
0.0000  0.2000  3.0000   % para 5  rrho_R2
3.0000  0.7500  2.0000   % para 6  kkappa
0.7000  0.1000  1.0000   % para 7  b
0.1500  0.0500  3.0000   % para 8  eeta_p
0.1500  0.0500  3.0000   % para 9  eeta_w
2.0000  0.7500  2.0000   % para 10 ttau
0.5000  0.2000  1.0000   % para 11 nnu_p
0.5000  0.2000  1.0000   % para 12 nnu_w
0.5000  0.2000  1.0000   % para 13 rrho_d
0.5000  0.2000  1.0000   % para 14 rrho_dL
0.5000  0.2000  1.0000   % para 15 rrho_A
0.5000  0.2000  1.0000   % para 16 rrho_mI
0.5000  0.2000  1.0000   % para 17 rrho_G
0.0010  2.0000  4.0000   % para 18 sigma_d
0.0010  2.0000  4.0000   % para 19 sigma_dL
0.0010  2.0000  4.0000   % para 20 sigma_A
0.0010  2.0000  4.0000   % para 21 sigma_mI
0.0010  2.0000  4.0000   % para 22 sigma_R
0.0010  2.0000  4.0000   % para 23 sigma_G
0.0010  2.0000  4.0000   % para 24 sigma_ppi (standard deviation of the measurement equation of inflation)
0.0010  2.0000  4.0000]; % para 25 sigma_ppi (standard deviation of the measurement equation of inflation)


% First Column: Prior mean
% Second Column: Prior Standard Deviation
% Third Column: Prior Shape:
%           1: BETA(mean,stdd)
%           2: GAMMA(mean,stdd)
%           3: NORMAL(mean,stdd)
%           4: INVGAMMA(s^2,nu)


end

