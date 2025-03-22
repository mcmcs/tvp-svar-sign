% function red_main 
% 1) Load data
% 2) Load prior 
% 3) Estimate split sample models 
% 4) Compute IRFS or BC decomposition 
close all; clear all; 

cucd=pwd; 
addpath([cucd,'\Chris Sims Gensys'],[cucd,'\Chris Sims csminwel']); 

% Model function 
funcmod =@modptjREDmain;    
% Minimization handle 
funcpost =@logpostSPLIT ;    
% Path to load data and prior 
loadpath=[cucd,'\baseline']; 

%% Load data 
cd(loadpath); 
[data,datanames]=xlsread('dataRED'); 
sample=data(:,1); 
Y     =data(:,2:end); 
datanames=datanames(2:end); 
clear data; 
cd(cucd); 

%% Initial parameter guess 
cd(loadpath); 
[param,parnames]=xlsread('initialguess'); 
cd(cucd); 

% Solve model with inial guess 
%[Gone,Cone,RRone,euone,SDXone,zmat,NY,NX,NETA,gev,ncof,Gtwo,RRtwo,eutwo,SDXtwo,Ctwo,Tbreak]=feval(funcmod,param); 

% Position entries calibrated
poscal=[ 2    13    31    33    36    37    38 ]'; 
% Position entries estimated 
posest=(1:length(param)); 
posest=find( ismember(posest(:),poscal(:)) == 0 ); 

% Transform initial guess to constrained minimization space 
xzero=boundsINV(param(posest)); 
logpostSPLIT(xzero,funcmod,param,posest,Y,1); 

%% Minimization using CSMINWEL 
[fh,xh,gh,H,itct,fcount,retcodeh] = csminwel(funcpost,xzero,.1*eye(length(xzero)),[],1e-6,1000,funcmod,param,posest,Y,1);

% Transform posterior from constrained minimization space to parameter
% space 
postmode=bounds(xh);

% Account for the transformation when computing Hessian 
JJ=jacobJPT(xh);
HH=JJ*H*JJ';
HH=0.5*(HH+HH'); 

try
    cd(loadpath); 
    save workspace; 
    xlswrite('tab_output',[parnames(posest) num2cell([postmode sqrt(diag(HH))])],'optimization'); 
    cd(cucd); 
catch 
    disp('Could not save'); 
end 
% ========================================================================
%% Begin Metropolis Algorithm 
% ========================================================================
% These are extremely conservative numbers for illustration purposes only
% Should use 4 chains with ntake=50,000 and nburn = 50,000. 
scaling=0.3; 
nchains=2;  
ntake  =2000; 
nburn  =2000;

drawsMat=zeros(length(xh),ntake+nburn,nchains); 

%% Generate starting values for each chain 
xStart      =zeros(length(xh),nchains); 
logPostStart=zeros(nchains   ,  1    ); 
Hscaled=scaling*scaling*HH; 
for ii=1:nchains
    logpostOLD=-1e10;
    while logpostOLD==-1e10;
        xguess=mvnrnd(postmode,4*Hscaled,1);
        logpostOLD=logpostSPLIT(xguess,funcmod,param,posest,Y,0); 
    end 
    drawsMat(:,1,ii)=xguess; 
    logPostStart(ii)=logpostOLD; 
end 
disp('Obtained starting values'); 

%% Begin each chain 
clockmat=zeros(nchains,1); 
acceptmat=zeros(nchains,1); 
for ii=1:nchains;     
    tic; 
    disp(sprintf('Begin chain %10.0f',ii));    
    count=0; 
    xOLD=drawsMat(:,1,ii); 
    logpostOLD=logPostStart(ii);     
    for jj=2:(ntake+nburn);
        if rem(jj,1000)==0;             
            disp(sprintf('Completed draws %8.0f with acceptance rate %2.4f',jj,count/jj)); 
            tiempo=toc; 
            clockmat(ii)=(tiempo/60) + clockmat(ii); 
            disp(sprintf('Run time for this block %3.2f',tiempo/60));             
            tic; 
        end         
        % Generate Candidate
        xcand=mvnrnd(drawsMat(:,jj-1,ii),Hscaled,1);
        logpostNEW=logpostSPLIT(xcand,funcmod,param,posest,Y,0);
        if logpostNEW > logpostOLD
            logpostOLD = logpostNEW;
            drawsMat(:,jj,ii)=xcand;
            count=count+1;
        else
            if rand(1)<exp(logpostNEW-logpostOLD);
                logpostOLD=logpostNEW;
                drawsMat(:,jj,ii)=xcand;
                count=count+1;
            else
                drawsMat(:,jj,ii)=drawsMat(:,jj-1,ii);
            end
        end                        
    end 
    acceptmat(ii)=count/(ntake+nburn);  % Acceptance Rate 
end 
disp(sprintf('Total acceptance rate %2.4f',mean(acceptmat))); 
%% 
try
    cd(loadpath); 
    save workspace; 
    cd(cucd); 
catch 
    disp('Could not save'); 
    
end 

%% Create table to summarize output 
tab_sum=zeros(length(xh),6); 
pvec=round([0.05 0.95]*ntake*nchains); 
for ii=1:length(xh); 
    xx=squeeze( drawsMat(ii,nburn+1:end,:) );
    xx=sort(xx(:));
    tab_sum(ii,1)=postmode(ii); 
    tab_sum(ii,2)=median(xx); 
    tab_sum(ii,3)=mean(xx); 
    tab_sum(ii,4)=std(xx); 
    tab_sum(ii,5:6)=xx(pvec); 
    clear xx; 
end   
cd(loadpath); 
xlswrite('tab_output',[parnames(posest) num2cell(tab_sum)],'mcmc'); 
cd(cucd); 
disp('End mcmc'); 