function simCluster(N,L,M,s2g,T,itcluster)

randn('seed',round(sum(1e5*clock)+itcluster));
rand('seed',round(sum(1e5*clock)+itcluster));

% T=200;
% %%Definimos el canal
% N=3;        % entradas
% L = 20;		% salidas
% s2g=1;
% M=1;
s2h=1;

generaMIMOrafagas;

BASEDIRini=['/export/clusterdata/franrruiz87/MIMO_N' num2str(N) '_L' num2str(L) '_M' num2str(M) '_s2g' num2str(10*s2g) '_T'  num2str(T)];
if ( isdir(BASEDIRini) == 0 )
    mkdir(BASEDIRini);
end

BASEDIR=['/export/clusterdata/franrruiz87/MIMO_N' num2str(N) '_L' num2str(L) '_M' num2str(M) '_s2g' num2str(10*s2g) '_T'  num2str(T) '/it'  num2str(itcluster)];
if ( isdir(BASEDIR) == 0 )
    mkdir(BASEDIR);
end

try
    load([BASEDIR '/SimFinalIt' num2str(itcluster) '.mat']);
catch err
    save([BASEDIR '/DatosIt' num2str(itcluster) '.mat']);
    s2x=2*s2g;
    alpha=.1;
    gamma1=.1;
    gamma2=10;
    beta=1;
    pii=0.2;
    lambda=1;
    D=L;


    Nsim=100;

    [Zest,Hest,MestIt,pX_ZIt,rest,nest] = MIMO_aprender(X,s2x,s2h,alpha,gamma1,gamma2,beta,pii,lambda,Nsim,itcluster,BASEDIR);

    save([BASEDIR '/SimFinalIt' num2str(itcluster) '.mat']);
end

