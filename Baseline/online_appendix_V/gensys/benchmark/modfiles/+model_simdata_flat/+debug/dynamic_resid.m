function [lhs, rhs] = dynamic_resid(y, x, params, steady_state)
T = NaN(0, 1);
lhs = NaN(10, 1);
rhs = NaN(10, 1);
lhs(1) = y(11);
rhs(1) = (-params(3))*(y(13)-y(22)-y(14))+y(21);
lhs(2) = y(12);
rhs(2) = y(22)*params(1)+y(11)*params(2)+y(17);
lhs(3) = y(13);
rhs(3) = params(4)*y(3)+(1-params(4))*(params(14)+y(12)*params(5)+y(11)*params(6))+y(15);
lhs(4) = y(14);
rhs(4) = params(14)+1/params(3)*params(7)*(params(9)-1)*y(16);
lhs(5) = y(15);
rhs(5) = params(8)*y(5)+params(11)*x(1);
lhs(6) = y(16);
rhs(6) = params(9)*y(6)+params(12)*x(2);
lhs(7) = y(17);
rhs(7) = params(10)*y(7)+params(13)*x(3);
lhs(8) = y(18);
rhs(8) = y(11)*100;
lhs(9) = y(19);
rhs(9) = 2+y(12)*400;
lhs(10) = y(20);
rhs(10) = y(13)*400;
end
