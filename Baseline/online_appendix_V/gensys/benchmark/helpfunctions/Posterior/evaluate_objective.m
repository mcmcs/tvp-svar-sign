function f = evaluate_objective(x)


global Y pmean pstdd pshape;

ttheta=x;%trans(x);

f = evaluate_posterior(ttheta,Y,pmean, pstdd, pshape)*-1;

end
