%% BCJR or FBs - T1

%% Uncomment one:
%carpeta = 'BCJR'; nombreFile = 'SimBCJRFinal_'; resultados = 'BCJR';
carpeta = 'FBsampling'; nombreFile = 'SimFBFinal_'; resultados = 'FBs';


Kvec=[2 3 4];
BER = zeros(11,3,100);
pOnOff = zeros(11,3,100);
n=0;
flag=0;
for N=[2:3 5]
    n=n+1;
    l=0;
    K=Kvec(n);
    for L=10
        for SNR=-5:5
             l=l+1;
            for T=400
                for itcluster=1:100
                    BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(flag) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];

                    load([BASEDIR '/' carpeta '/' nombreFile num2str(itcluster) '.mat']);

                    simb=simb(:,K:end);
                    idxSimb = find(simb~=0);

                    BER(l,n,itcluster) =sum(sum(Zest~=simb))/(N*T);%sum(Zest(idxSimb)~=simb(idxSimb))/length(idxSimb);
                    pOnOff(l,n,itcluster) = sum(sum(abs(Zest)~=abs(simb)))/(N*T);
                end                            
            end
        end
    end
end


save (['resultados' resultados '_SNR'], 'BER','pOnOff');
