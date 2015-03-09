clear all;
N=2;
n=1;
if N==2
load('resultadosBCJR_SNR.mat');
BERBCJR=BER;
ADERBCJR=pOnOff;
end
load('resultadosFBs_SNR.mat');
BERFB=BER;
ADERFB=pOnOff;

load('resultadosMIMO_SNR.mat');


SNRvec=-5:5;
ber=zeros(length(SNRvec),1);
mse=zeros(length(SNRvec),1);
ader=zeros(length(SNRvec),1);
dep=zeros(length(SNRvec),1);
udep=zeros(length(SNRvec),1);

l=0;
for L=SNRvec
    l=l+1;
    idxA=1:100;%find(Nfin(l,n,:)==N);%find(Nfin(l,n,:)==N & Mfin(l,n,:)==1);
    ber(l)=sum(BER(l,n,idxA)-BERFB(l,n,idxA))/length(idxA);
    e(l)=std(BER(l,n,idxA)-BERFB(l,n,idxA));
    
    mse(l)=sum(MSE(l,n,idxA))/(length(idxA)*L);
    ader(l)=sum(pOnOff(l,n,idxA)-ADERFB(l,n,idxA))/length(idxA);
    eAder(l)=std(pOnOff(l,n,idxA)-ADERFB(l,n,idxA));
    dep(l)=1-sum(Mfin(l,n,:) & Nfin(l,n,:)==N)/100;
    udep(l)=1-sum(Nfin(l,n,:)==N)/100;
    
    if N==2
        ber2(l)=sum(BER(l,n,idxA)-BERBCJR(l,n,idxA))/length(idxA);
        e2(l)=std(BER(l,n,idxA)-BERBCJR(l,n,idxA));
        ader2(l)=sum(pOnOff(l,n,idxA)-ADERBCJR(l,n,idxA))/length(idxA);
        eAder2(l)=std(pOnOff(l,n,idxA)-ADERBCJR(l,n,idxA));
    end
    
    
end

errorbar(SNRvec,ber,e,'r');
if N==2
hold on, errorbar(SNRvec,ber2,e2);
end
figure,errorbar(SNRvec,ader,eAder,'r');
if N==2
hold on, errorbar(SNRvec,ader2,eAder2);
end
