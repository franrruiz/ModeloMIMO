function simFer(itH,itcluster,SNR)
BASEDIR=['/export/clusterdata/ivalera/modMIMOFerN3T400_L10/SNR' num2str(SNR) '_itH'  num2str(itH)];
%BASEDIR=['./flag' num2str(flag) '_flagChi' num2str(flagChi) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
if ( isdir(BASEDIR) == 0 )
    mkdir(BASEDIR);
end

try
    load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
catch err    

    randn('seed',itcluster+round(sum(1e5*clock)));
    rand('seed',itcluster+round(sum(1e5*clock)));


    %%Definimos el canal
    T=400;
    L=10;
    N=3;        % entradas
    Mvec=[3 2 1];
    M=max(Mvec);
    s2g= 10^(-SNR/10);
   
    chi=1;%1/s2g;


    load('symbols.mat');
    load(['H_it' num2str(itH) '.mat']);
    
    X=zeros(L,T);

    for t=1:T
       y=simb(:,t:t+M-1); 
       X(:,t)=H*y(:)+sqrt(s2g)*randn(L,1);
    end

    tau=1;
    nu=.1;
    alpha1=.1;
    alpha2=1;
    gamma1=.1;
    gamma2=10;
    beta=1;
    pii=0.2;
    lambda=1;
    D=L;

    Nsim=1000;

    [Zest,Hest,MestIt,pX_ZIt,chiIt,rest,nest] = MIMO_aprender_H_s2x(X,chi,tau,nu,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, 0);

    save([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
end
% 
% exit
