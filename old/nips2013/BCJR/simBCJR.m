function simBCJR(L,SNR,T,itcluster)
BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN2T' num2str(T) '/flag' num2str(0) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
%BASEDIR=['./flag' num2str(flag) '_flagChi' num2str(flagChi) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
if ( isdir([BASEDIR '/BCJR']) == 0 )
    mkdir([BASEDIR '/BCJR']);
end

try
    load([BASEDIR '/BCJR/SimBCJRFinal_' num2str(itcluster) '.mat']);
catch err    
    
    load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
    randn('seed',itcluster+round(sum(1e5*clock)));
    rand('seed',itcluster+round(sum(1e5*clock)));

    
    Zest=bcjr(X,H,s2g);

    save([BASEDIR '/BCJR/SimBCJRFinal_' num2str(itcluster) '.mat'],'Zest', 'simb');
end

