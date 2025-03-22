function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = model_simdata_flat.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(10, 1);
    residual(1) = (y(11)) - ((-params(3))*(y(13)-y(22)-y(14))+y(21));
    residual(2) = (y(12)) - (y(22)*params(1)+y(11)*params(2)+y(17));
    residual(3) = (y(13)) - (params(4)*y(3)+(1-params(4))*(params(14)+y(12)*params(5)+y(11)*params(6))+y(15));
    residual(4) = (y(14)) - (params(14)+1/params(3)*params(7)*(params(9)-1)*y(16));
    residual(5) = (y(15)) - (params(8)*y(5)+params(11)*x(1));
    residual(6) = (y(16)) - (params(9)*y(6)+params(12)*x(2));
    residual(7) = (y(17)) - (params(10)*y(7)+params(13)*x(3));
    residual(8) = (y(18)) - (y(11)*100);
    residual(9) = (y(19)) - (2+y(12)*400);
    residual(10) = (y(20)) - (y(13)*400);
end
