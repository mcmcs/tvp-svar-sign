function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(15)=params(8)*y(5)+params(11)*x(1);
  y(16)=params(9)*y(6)+params(12)*x(2);
  y(17)=params(10)*y(7)+params(13)*x(3);
  y(14)=params(14)+1/params(3)*params(7)*(params(9)-1)*y(16);
end
