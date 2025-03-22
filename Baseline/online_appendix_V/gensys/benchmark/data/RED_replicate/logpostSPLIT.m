function logpost=logpostSPLIT(pzero,funcmod,param,posest,Y,flag_min)
% =========================================================================
% Date: February 22 2010
% This function computes the value of the posterior density for the model
% in Justiniano, Primiceri and Tambalotti (RED 2010), allowing for split
% samples 
% 
% FUNCMOD   Function handle 
% PARAM     Full parameter vector 
% POSEST    Entries in parameter vector that are estimated 
% PZERO     Estimated parameters 
% Y         [T,n] data matrix 
% Set FLAG_MIN == 1 to use in Minimization 
%     FLAG_MIN == 0 to use in Maximization (e.g. MCMC) 
% =========================================================================

%% If minimization, transform the unconstrained into a constrained minimization
if flag_min==1
    param(posest)=bounds(pzero);
    logpost=1e10;
else
    param(posest)=pzero;
    logpost= -1e10;
end

%% Compute prior density, requires full parameter vector
[logprior,flag_ok]=logpriorJPT(param);
if flag_ok==0;return;end  

%% Solve model 
[G1,CC1,RR,eu,SDX,ZZ,~,~,~,~,~,G2,RR2,eu2,SDX2,CC2,Tbreak]=feval(funcmod,param);
if isequal(eu,[1;1])==0 || isequal(eu2,[1;1])==0
     return;
end

[T,nn]=size(Y);

%% Model-based demeaning 
ss=size(G1,1);
if any(CC1~=0)==1
    Y(1:Tbreak,:)=Y(1:Tbreak,:)-repmat((ZZ*CC1)',[Tbreak 1]);
end 
if any(CC2~=0)==1;
    Y(Tbreak+1:end,:)=Y(Tbreak+1:end,:)-repmat((ZZ*CC2)',[(T-Tbreak) 1]);
end 
Y=Y'; 

%% Initialize Kalman filter 
MM=RR*(SDX');
try
    pshat=disclyap_fast(G1, MM*(MM'));
catch
    disp('Warning! Problems with initialization in GENSYSPOST');
    return;
end
shat=zeros(ss,1);lht=zeros(T,1);

%% Kalman filter Loop 
% Sample 1
for ii=1:Tbreak;
    [shat,pshat,lht(ii)]=kf(Y(:,ii),ZZ,shat,pshat,G1,MM);
end
% Sample 2 
MM=RR2*(SDX2');
for ii=(Tbreak+1):T;
    [shat,pshat,lht(ii)]=kf(Y(:,ii),ZZ,shat,pshat,G2,MM);
end
logpost=-((T*nn*0.5)*(log(2*pi)))+(sum(lht)+logprior);
if flag_min==1
    logpost=-1*logpost;
end 