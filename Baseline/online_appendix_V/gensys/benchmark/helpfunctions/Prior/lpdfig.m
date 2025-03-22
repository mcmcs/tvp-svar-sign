
function [ldens] = lpdfig(x,a,b)
% log INVERSE GAMMA
ldens = - lngam(a) + a*log(b) - (a+1)*log(x) - b/x;
end

% Zellner (1971) page 370. Equation A.37a

