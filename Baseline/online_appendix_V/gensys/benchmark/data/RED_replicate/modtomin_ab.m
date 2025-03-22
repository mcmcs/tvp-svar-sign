function theta=modtomin_ab(x,a,b) 
% Transform a parameter in the [a,b] interval to (-inf,inf) 
% Inverse function is mintomod_ab.m 
c=2*( x-0.5*(a+b) )/(b-a); 
theta=c/sqrt(1-c*c);