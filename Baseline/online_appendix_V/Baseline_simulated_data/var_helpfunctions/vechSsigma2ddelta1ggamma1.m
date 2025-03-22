function y = vechSsigma2ddelta1ggamma1(x,info)


Ssigma = inv_vec(x,info.nvar);

D          = sqrt(diag(diag(Ssigma)));
matrixC    = (D\Ssigma)/D;

logmatrixC = logm(matrixC);
ggamma     = logmatrixC(logical(tril(logmatrixC,-1)));

y          = [2*log(diag(D));ggamma;];
end

