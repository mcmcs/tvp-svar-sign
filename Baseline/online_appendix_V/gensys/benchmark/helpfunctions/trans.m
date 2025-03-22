function [ttheta] = trans(ttheta)

% this function transforms variables from max to model


transspec = get_transspec();


npara = size(ttheta,1);

i = 1;
while i <= npara;

   a = transspec(i,2);
   b = transspec(i,3);
   c = transspec(i,4);

   if transspec(i,1) == 1;
     
     ttheta(i) = (a+b)/2 + 0.5*(b-a)*c*ttheta(i)/sqrt(1+c^2*ttheta(i)^2);
     

   elseif transspec(i,1) ==2;
     ttheta(i) = a + exp(c*(ttheta(i)-b));
   end;

   i = i+1;

end;


