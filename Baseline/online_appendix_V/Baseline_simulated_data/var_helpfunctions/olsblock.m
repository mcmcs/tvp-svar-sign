% Obtained from Del Negro and Primiceri's (2010)
% FORECASTING by OLS
function r=olsblock(y,x) 

if (nargin ~= 2); error('Wrong # of arguments to ols');
else
 [nobs nvar] = size(x); [nobsy junk] = size(y);
 if (nobs ~= nobsy); error('x and y must have same # obs in ols');
 end;
end;
r.nobs=nobs;
r.nvar=nvar;

r.bhatols=(x'*x)\(x'*y);
r.yhatols=x*r.bhatols;
r.resols=y-r.yhatols;
r.sigmahatols = (r.resols'*r.resols)/(nobs-nvar);

r.sig2hatols=sum(r.resols.^2)/(nobs-nvar);
%r.sigbhatols=r.sig2hatols*(inv(x'*x));
r.XX=(x'*x);
%r.R2=var(r.yhatols)/var(y);


