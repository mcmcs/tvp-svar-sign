function J=jacobJPT(param);
% Entry vector is only for estimated parameters 
J=zeros(length(param),1);
J(1)=bound01prime(param(1));     
J(2)=bound01prime(param(2));     
J(3)=bound01prime(param(3));     
J(4)=1;
J(5)=1; 
J(6)=bound01prime(param(6));  
J(7)=1; 
J(8)=1; 
J(9)=1;       
J(10)=1;   
J(11)=bound0prime(param(11));
J(12)=bound0prime(param(12));
J(13)=bound01prime(param(13));
J(14)=bound01prime(param(14));
J(15)=bound0prime(param(15));
J(16)=bound0prime(param(16));
J(17)=1; 
J(18)=1;       
J(19)=bound0Bprime(param(19),0,.99);
J(20)=bound0Bprime(param(20),0,.99);
J(21)=bound0Bprime(param(21),0,.99);
J(22)=bound0Bprime(param(22),0,.99);
J(23)=bound0Bprime(param(23),0,.99);
J(24)=bound0Bprime(param(24),0,.99);
J(25)=bound0Bprime(param(25),0,.99);
J(26)=bound0Bprime(param(26),0,.99);
J(27)=bound0Bprime(param(27),0,.99);
J(28)=1;
J(29)=bound0Bprime(param(29),0,0.99);
J(30)=1; 
J(31)=1; 
J(32)=bound0prime(param(32));
J(33)=bound0prime(param(33)); 
J(34)=bound0prime(param(34)); 
J(35)=bound0prime(param(35));
J(36)=bound0prime(param(36)); 
J(37)=bound0prime(param(37)); 
J(38)=bound0prime(param(38));
J(39)=bound0prime(param(39)); 

if any(J==0)
    error('Zero entries on diagonal') 
end
if any(imag(J)~=0) 
    error('Complex entries on diagonal'); 
end 
if any(isinf(abs(J))~=0) 
    error('Inf entries on diagonal') 
end 
J=diag(J); 

function y = bound01prime(x);
y = exp(x)/(1+exp(x))^2;

function y=bound0Bprime(x,a,b);
y=0.5*(b-a)/((1+x*x)^1.5);

function y = bound0prime(x);
y=exp(x);