function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model_simdata.static_resid_tt(T, y, x, params);
end
residual = zeros(10, 1);
    residual(1) = (y(1)) - (y(1)+(-params(3))*(y(3)-y(2)-y(4)));
    residual(2) = (y(2)) - (y(2)*params(1)+y(1)*params(2)+y(7));
    residual(3) = (y(3)) - (y(3)*params(4)+(1-params(4))*(params(14)+y(2)*params(5)+y(1)*params(6))+y(5));
    residual(4) = (y(4)) - (params(14)+1/params(3)*params(7)*(params(9)-1)*y(6));
    residual(5) = (y(5)) - (y(5)*params(8)+params(11)*x(1));
    residual(6) = (y(6)) - (params(9)*y(6)+params(12)*x(2));
    residual(7) = (y(7)) - (y(7)*params(10)+params(13)*x(3));
    residual(8) = (y(8)) - (y(1)*100);
    residual(9) = (y(9)) - (2+y(2)*400);
    residual(10) = (y(10)) - (y(3)*400);

end
