clear all;
N=5;
n=3;
load('resultadosMIMO_L.mat');


Lvec=5:3:20;
ber=zeros(length(Lvec),1);
mse=zeros(length(Lvec),1);
ader=zeros(length(Lvec),1);
dep=zeros(length(Lvec),1);
udep=zeros(length(Lvec),1);

l=0;
for L=Lvec
    l=l+1;
    idxA=find(Nfin(l,n,:)==N & Mfin(l,n,:)==1);
    ber(l)=sum(BER(l,n,idxA))/length(idxA);
    upp(l)=std(BER(l,n,idxA));
    low(l)=min(upp(l),ber(l));
    mse(l)=sum(MSE(l,n,idxA))/(length(idxA)*L);
    ader(l)=sum(pOnOff(l,n,idxA))/length(idxA);
    dep(l)=1-sum(Mfin(l,n,:) & Nfin(l,n,:)==N)/100;
    udep(l)=1-sum(Nfin(l,n,:)==N)/100;
end
load('resultadosFBs_L.mat');
berFB=zeros(length(Lvec),1);
aderFB=zeros(length(Lvec),1);

l=0;
for L=Lvec
    l=l+1;
    berFB(l)=sum(BER(l,n,:))/100;  
    uppFB(l)=std(BER(l,n,:));
    lowFB(l)=min(uppFB(l),berFB(l));
    aderFB(l)=sum(pOnOff(l,n,:))/100;
end

if N==2
    load('resultadosBCJR_L.mat');
    berBCJR=zeros(length(Lvec),1);
    aderBCJR=zeros(length(Lvec),1);

    l=0;
    for L=Lvec
        l=l+1;
        berBCJR(l)=sum(BER(l,n,:))/100;  
        uppBCJR(l)=std(BER(l,n,:));
        lowBCJR(l)=min(uppBCJR(l),berBCJR(l));
        aderBCJR(l)=sum(pOnOff(l,n,:))/100;
    end
    figure, errorbar(Lvec,berBCJR,lowBCJR, uppBCJR);
    hold on, errorbar(Lvec,ber,low,upp,'r');
    errorbar(Lvec,berFB,lowFB,uppFB,'k');

    %figure,plot(Lvec, ber, Lvec, berFB, Lvec, berBCJR);
    %figure,plot(Lvec, ader, Lvec, aderFB, Lvec, aderBCJR);
else
    %figure,plot(Lvec, ber, Lvec, berFB);
    %figure,plot(Lvec, ader, Lvec, aderFB);
    errorbar(Lvec,ber,low,upp,'r');
    hold on, errorbar(Lvec,berFB,lowFB,uppFB,'k');
    
    
end
figure,plot(Lvec, dep);
figure,plot(Lvec, udep);
figure,plot(Lvec, mse);