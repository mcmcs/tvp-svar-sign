function [ttheta] = invtrans(ttheta)

% this procedure transforms variables from model to max
% Note that the max parameter do not contain para 8 and 9


transspec = get_transspec();

npara = size(ttheta,1);

i = 1;
while i <= npara;

   a = transspec(i,2);
   b = transspec(i,3);
   c = transspec(i,4);

   if transspec(i,1) == 1;
     cx = 2*(ttheta(i)-(a+b)/2)/(b-a);
     ttheta(i) = (1/c)*cx/sqrt(1-cx^2);
   elseif transspec(i,1) == 2;
     ttheta(i) = b + (1/c)*log(ttheta(i)-a);
   end;

   i = i+1;


end



