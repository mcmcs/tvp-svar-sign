function [y, T] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(10)=y(3)*400;
  y(9)=2+y(2)*400;
  y(8)=y(1)*100;
end
