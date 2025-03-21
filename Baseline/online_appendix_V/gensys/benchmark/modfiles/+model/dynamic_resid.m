function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(7, 1);
lhs = y(5);
rhs = (-params(3))*(y(7)-y(13)-y(8))+y(12);
residual(1) = lhs - rhs;
lhs = y(6);
rhs = y(13)*params(1)+y(5)*params(2)+y(11);
residual(2) = lhs - rhs;
lhs = y(7);
rhs = params(4)*y(1)+(1-params(4))*(params(14)+y(6)*params(5)+y(5)*params(6))+y(9);
residual(3) = lhs - rhs;
lhs = y(8);
rhs = params(14)+1/params(3)*params(7)*(params(9)-1)*y(10);
residual(4) = lhs - rhs;
lhs = y(9);
rhs = params(8)*y(2)+params(11)*x(it_, 1);
residual(5) = lhs - rhs;
lhs = y(10);
rhs = params(9)*y(3)+params(12)*x(it_, 2);
residual(6) = lhs - rhs;
lhs = y(11);
rhs = params(10)*y(4)+params(13)*x(it_, 3);
residual(7) = lhs - rhs;

end
