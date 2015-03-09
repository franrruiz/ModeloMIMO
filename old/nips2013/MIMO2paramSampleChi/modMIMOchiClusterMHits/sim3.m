function sim2(L,SNR,T,itcluster,flag)
BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN5T400/flag' num2str(flag) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
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
    N=5;        % entradas
    Mvec=[4 4 3 2 2];
    s2h=1;
    s2g= 10^(-SNR/10);
   
    chi=1;%1/s2g;


    generaMIMOrafagas2;

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

    [Zest,Hest,MestIt,pX_ZIt,chiIt,rest,nest] = MIMO_aprender_H_s2x(X,chi,tau,nu,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, flag);

    save([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
end
% 
% exit
