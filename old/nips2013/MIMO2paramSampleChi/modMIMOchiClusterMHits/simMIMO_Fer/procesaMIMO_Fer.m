for itH=1:20
    for SNR=-5:5
        procesaMIMO_Fer_cluster(itH,SNR);
    end
end


%% Pocresa SNR
Nfin = zeros(11,20,10);
Mfin = zeros(11,20,10);
MSE = zeros(11,20,10);
BER = zeros(11,20,10);
MEM = zeros(11,20,10,3);
pOnOff = zeros(11,20,10);
flag=0;
for itH=1:20
    n=itH;
    for L=10
        l=0;
        for SNR=-5:5
            l=l+1;
            for T=400
                BASEDIR=['/export/clusterdata/ivalera/modMIMOFerN3T400_L10/SNR' num2str(SNR) '_itH'  num2str(itH)];

                load([BASEDIR '/../Resultados/SNR' num2str(SNR) '_itH'  num2str(itH) '.mat']);


                Mfin(l,n,:) = Mfin_;
                Nfin(l,n,:) = Nfin_;

                BER(l,n,:) = BER_;
                pOnOff(l,n,:) = pOnOff_;
                MSE(l,n,:) = MSE_;
                
                MEM(l,n,:,:) = mem_;
            end
        end

    end
end



save ('resultadosMIMO_Fer.mat', 'BER','MSE','Mfin','Nfin','pOnOff','MEM')
