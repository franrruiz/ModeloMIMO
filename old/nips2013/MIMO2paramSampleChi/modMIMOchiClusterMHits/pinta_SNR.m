close all;clear all;
N=2;
n=1;
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
    idxA=find(Nfin(l,n,:)==N & Mfin(l,n,:)==1);
    %idxA=1:100;
    ber(l)=sum(BER(l,n,idxA))/length(idxA);
    upp(l)=std(BER(l,n,idxA));
    low(l)=min(upp(l),ber(l));
    mse(l)=sum(MSE(l,n,idxA))/(length(idxA)*10);
    ader(l)=sum(pOnOff(l,n,idxA))/length(idxA);
    dep(l)=1-sum(Mfin(l,n,:) & Nfin(l,n,:)==N)/100;
    udep(l)=1-sum(Nfin(l,n,:)==N)/100;
end

load('resultadosFBs_SNR.mat');
berFB=zeros(length(SNRvec),1);
aderFB=zeros(length(SNRvec),1);

l=0;
for L=SNRvec
    l=l+1;
    berFB(l)=sum(BER(l,n,idxA))/length(idxA);  
    uppFB(l)=std(BER(l,n,idxA));
    lowFB(l)=min(uppFB(l),berFB(l));
    aderFB(l)=sum(pOnOff(l,n,idxA))/length(idxA);
end

if N==2
    load('resultadosBCJR_SNR.mat');
    berBCJR=zeros(length(SNRvec),1);
    aderBCJR=zeros(length(SNRvec),1);

    l=0;
    for L=SNRvec
        l=l+1;
        berBCJR(l)=sum(BER(l,n,:))/100;  
        uppBCJR(l)=std(BER(l,n,:));
        lowBCJR(l)=min(uppBCJR(l),berBCJR(l));
        aderBCJR(l)=sum(pOnOff(l,n,:))/100;
    end
    
    figure,errorbar(SNRvec,ber,low,upp,'r');
    hold on, errorbar(SNRvec,berFB,lowFB,uppFB,'k');
    errorbar(SNRvec,berBCJR,lowBCJR, uppBCJR);
    legend('IFU-FSM', 'FW-BW sampling','BCJR'), ylabel('BER'),xlabel('-10 log(\sigma_n^2)')
%     figure,errorbar(SNRvec,ber,low,upp,'r');
%     hold on, errorbar(SNRvec,berFB,lowFB,uppFB,'k');
%     errorbar(SNRvec,berBCJR,lowBCJR, uppBCJR);
%     legend('IFU-FSM', 'FW-BW sampling','BCJR'), ylabel('ADER'),xlabel('-10 log(\sigma_n^2)')
else
    figure,errorbar(SNRvec,ber,low,upp,'r');
    hold on, errorbar(SNRvec,berFB,lowFB,uppFB,'k');
    legend('IFU-FSM', 'FW-BW sampling'), ylabel('BER'),xlabel('-10 log(\sigma_n^2)')
    
end
figure,plot(SNRvec, dep);
legend('IFU-FSM'), ylabel('DEP'),xlabel('-10 log(\sigma_n^2)')
figure,plot(SNRvec, udep);
legend('IFU-FSM'), ylabel('UDEP'),xlabel('-10 log(\sigma_n^2)')
figure,plot(SNRvec, mse);
legend('IFU-FSM'), ylabel('MSE'),xlabel('-10 log(\sigma_n^2)')
