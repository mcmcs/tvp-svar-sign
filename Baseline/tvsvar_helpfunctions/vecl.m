function y = vecl(x)
% extract lower off-diagonal elements of x
% y = x(find(tril(true(size(x)), -1)));
y = x( tril(true(size(x)), -1) == 1);