function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = model_simdata_flat.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(10, 1);
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
