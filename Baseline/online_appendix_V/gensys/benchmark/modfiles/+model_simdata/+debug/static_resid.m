function [lhs, rhs] = static_resid(y, x, params)
T = NaN(0, 1);
lhs = NaN(10, 1);
rhs = NaN(10, 1);
lhs(1) = y(1);
rhs(1) = y(1)+(-params(3))*(y(3)-y(2)-y(4));
lhs(2) = y(2);
rhs(2) = y(2)*params(1)+y(1)*params(2)+y(7);
lhs(3) = y(3);
rhs(3) = y(3)*params(4)+(1-params(4))*(params(14)+y(2)*params(5)+y(1)*params(6))+y(5);
lhs(4) = y(4);
rhs(4) = params(14)+1/params(3)*params(7)*(params(9)-1)*y(6);
lhs(5) = y(5);
rhs(5) = y(5)*params(8)+params(11)*x(1);
lhs(6) = y(6);
rhs(6) = params(9)*y(6)+params(12)*x(2);
lhs(7) = y(7);
rhs(7) = y(7)*params(10)+params(13)*x(3);
lhs(8) = y(8);
rhs(8) = y(1)*100;
lhs(9) = y(9);
rhs(9) = 2+y(2)*400;
lhs(10) = y(10);
rhs(10) = y(3)*400;
end
