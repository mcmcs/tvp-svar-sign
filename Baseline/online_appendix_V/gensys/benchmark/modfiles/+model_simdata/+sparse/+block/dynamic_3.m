function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(20)=y(13)*400;
  y(19)=2+y(12)*400;
  y(18)=y(11)*100;
end
