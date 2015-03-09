function simMIMO_FB(L,SNR,T,itcluster,N)
BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(0) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
%BASEDIR=['./flag' num2str(flag) '_flagChi' num2str(flagChi) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
if ( isdir([BASEDIR '/MIMOFB']) == 0 )
    mkdir([BASEDIR '/MIMOFB']);
end

try
    load([BASEDIR '/MIMOFB/SimFinal_' num2str(itcluster) '.mat']);
catch err    
    
    
    load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
    randn('seed',itcluster+round(sum(1e5*clock)));
    rand('seed',itcluster+round(sum(1e5*clock)));

    Nsim=10000;
    alpha=1;
    beta1=1;
    chi=1;
    try
        [Zest,Hest,MestIt,pX_ZIt,chiIt,rest,nest] = MIMO_aprender_H_s2x(X,chi,tau,nu,alpha,gamma1,gamma2,beta1,Nsim);
        save([BASEDIR '/MIMOFB/SimFinal_' num2str(itcluster) '.mat']);
    catch err2
        disp(err2);
    end
        

    
end



