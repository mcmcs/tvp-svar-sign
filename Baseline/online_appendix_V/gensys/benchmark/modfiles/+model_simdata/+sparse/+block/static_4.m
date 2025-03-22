function [y, T] = static_4(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(4)=params(14)+1/params(3)*params(7)*(params(9)-1)*y(6);
end
