function simFBsampling(L,SNR,T,itcluster,N)
BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(0) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
%BASEDIR=['./flag' num2str(flag) '_flagChi' num2str(flagChi) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
if ( isdir([BASEDIR '/FBsampling']) == 0 )
    mkdir([BASEDIR '/FBsampling']);
end

try
    load([BASEDIR '/FBsampling/SimFBFinal_' num2str(itcluster) '.mat']);
catch err    
    
    load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
    randn('seed',itcluster+round(sum(1e5*clock)));
    rand('seed',itcluster+round(sum(1e5*clock)));

    Nsim=10000;
    
    Zest=FBsampling(X,Mvec,H,Nsim,s2g);

    save([BASEDIR '/FBsampling/SimFBFinal_' num2str(itcluster) '.mat'],'Zest', 'simb');
end



