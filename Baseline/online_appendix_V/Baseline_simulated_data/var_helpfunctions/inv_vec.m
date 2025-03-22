function y = inv_vec(x,n)

Dn=duplication_matrix(n);

y=kron(vec(eye(n))',eye(n))*kron(eye(n),vec(Dn*x));