%% BCJR or FBs - T1

%% Uncomment one:
carpeta = 'BCJR'; nombreFile = 'SimBCJRFinal_'; resultados = 'BCJR';
%carpeta = 'FBsampling'; nombreFile = 'SimFBFinal_'; resultados = 'FBs';


Kvec=[2 3 4];
BER = zeros(6,3,100);
pOnOff = zeros(6,3,100);
n=0;
flag=0;
for N=2%[2:3 5]
    n=n+1;
    l=0;
    K=Kvec(n);
    for L=5:3:20
       l=l+1;
        for SNR=0
            for T=400
                for itcluster=1:100
                    BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(flag) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];

                    load([BASEDIR '/' carpeta '/' nombreFile num2str(itcluster) '.mat']);

                    simb=simb(:,K:end);
                    idxSimb = find(simb~=0);

                    BER(l,n,itcluster) =sum(Zest(:)~=simb(:))/(N*T);
                    pOnOff(l,n,itcluster) = sum(sum(abs(Zest)~=abs(simb)))/(N*T);
                end                            
            end
        end
    end
end


save (['resultados' resultados '_L'], 'BER','pOnOff');
