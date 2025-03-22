function [ldens] = lpdfgam(x,a,b)
% log GAMMA PDF
ldens = -lngam(a) -a*log(b)+ (a-1)*log(x) -x/b ;
end
