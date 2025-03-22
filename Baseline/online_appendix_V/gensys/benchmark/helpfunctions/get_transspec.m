function  transspec = get_transspec()

% column 1: Type of transformation  (do not change this)
% column 2: Lower bound of estimated parameters
% column 3: Upper bound of estimated parameters
% column 4: Parameter associated with the transformation (do not change this)


transspec=[
    2   1E-5     0       1 % para 1  ppsi_ppi response to inflation
    2   1E-5     0       1 % para 2  ppsi_gy  response to output growth
    2   1E-5     0       1 % para 3  ppsi_y
    1   -2       2       1 % para 4  rrho_R1
    1   -2       2       1 % para 5  rrho_R2
    2   1E-5     0       1 % para 6  kkappa
    1   1E-5     0.9999  1 % para 7  b
    1   -2       2       1 % para 8  eeta_p
    1   -2       2       1 % para 9  eeta_w
    2   1E-5     0       1 % para 10 ttau
    1   1E-5     0.9999  1 % para 12 nnu_p
    1   1E-5     0.9999  1 % para 13 nnu_w
    1   1E-5     0.9999  1 % para 14 rrho_d
    1   1E-5     0.9999  1 % para 15 rrho_dL
    1   1E-5     0.9999  1 % para 16 rrho_A
    1   1E-5     0.9999  1 % para 17 rrho_mI
    1   1E-5     0.9999  1 % para 19 rrho_G
    2   1E-5     0       1 % para 20 std shock d
    2   1E-5     0       1 % para 21 std shock d_L
    2   1E-5     0       1 % para 22 std shock A
    2   1E-5     0       1 % para 23 std shock m_I
    2   1E-5     0       1 % para 25 std shock R
    2   1E-5     0       1 % para 26 std shock G
    2   1E-5     0       1 % para 27 std shock measurement equation for inflation
    2   1E-5     0       1 % para 28 std shock measurement equation for wage growth
    ];

end