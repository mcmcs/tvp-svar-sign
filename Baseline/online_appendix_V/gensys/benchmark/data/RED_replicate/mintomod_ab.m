function x=mintomod_ab(theta,a,b) 
% Transform a parameter in the (-inf,inf) interval to [a,b]
% This is the inverse function to modtomin_ab.m 
x=0.5*(a+b)+0.5*(b-a)*theta/sqrt(1+theta*theta);