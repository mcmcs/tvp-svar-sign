

function [lng] = lngam(x)
% This procedure computes the log-gamma function


  if x <= 3
     lng= log(gamma(x));
  else
     k = floor(x-2);
     lng= sum(log(seqa(x-1,-1,k))) + log(gamma(x-k));
  end

end
