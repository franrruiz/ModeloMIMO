%% Pocresa L

Nfin = zeros(6,3,100);
Mfin = zeros(6,3,100);
MSE = zeros(6,3,100);
BER = zeros(6,3,100);
MEM = zeros(6,3,100,5);
pOnOff = zeros(6,3,100);
n=0;
flag=0;
for N=[2:3 5]
    n=n+1;
    l=0;
    for L=5:3:20
       l=l+1;
        for SNR=0
            for T=400
                BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(flag) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
%                     if ( isdir(BASEDIR) == 0 )
%                         mkdir(BASEDIR);
%                     end

                load([BASEDIR '/../ResultadosBER/N' num2str(N) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T) '.mat']);


                Mfin(l,n,:) = Mfin_;
                Nfin(l,n,:) = Nfin_;

                BER(l,n,:) = BER_;
                pOnOff(l,n,:) = pOnOff_;
                MSE(l,n,:) = MSE_;
                MEM(l,n,:,1:N ) = mem_;
            end
        end

    end
end



save ('resultadosMIMO_L.mat', 'BER','MSE','Mfin','Nfin','pOnOff','MEM')
