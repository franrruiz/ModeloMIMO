close all;clear all;
N=3;
load('resultadosMIMO_Fer.mat');


SNRvec=-5:5;
ber=zeros(length(SNRvec),1);
mse=zeros(length(SNRvec),1);
ader=zeros(length(SNRvec),1);
dep=zeros(length(SNRvec),1);
udep=zeros(length(SNRvec),1);
mem1=zeros(length(SNRvec),1);
mem2=zeros(length(SNRvec),1);
mem3=zeros(length(SNRvec),1);

l=0;
sum_idxA=0;
%itH=2;
%n=itH;
vec_itH=1:20;
for L=SNRvec
    l=l+1;
    for itH=vec_itH
        n=itH;
        %idxA=find(Nfin(l,n,:)==N);
        idxA=find(Nfin(l,n,:)==N & Mfin(l,n,:)==1);
        %idxA=1:10;
        sum_idxA=sum_idxA+length(idxA);
        ber(l)=ber(l)+sum(BER(l,n,idxA));
        %upp(l)=std(BER(l,n,idxA));
        %low(l)=min(upp(l),ber(l));
        mse(l)=mse(l)+sum(MSE(l,n,idxA))/10;
        ader(l)=ader(l)+sum(pOnOff(l,n,idxA));
        dep(l)=dep(l)+1-sum(Mfin(l,n,:) & Nfin(l,n,:)==N)/10;
        udep(l)=udep(l)+1-sum(Nfin(l,n,:)==N)/10;
        
        idxA=find(Nfin(l,n,:)==N);
        mem1(l)=mem1(l)+sum(MEM(l,n,idxA,1));
        mem2(l)=mem2(l)+sum(MEM(l,n,idxA,2));
        mem3(l)=mem3(l)+sum(MEM(l,n,idxA,3));
    end
end
dep=dep/length(vec_itH);
udep=udep/length(vec_itH);
ber=ber/sum_idxA;
ader=ader/sum_idxA;
mse=mse/sum_idxA;
mem1=mem1./sum(sum(Nfin(:,vec_itH,:)==3,2),3);
mem2=mem2./sum(sum(Nfin(:,vec_itH,:)==3,2),3);
mem3=mem3./sum(sum(Nfin(:,vec_itH,:)==3,2),3);


%load('resultadosFBs_SNR.mat');
%berFB=zeros(length(SNRvec),1);
%aderFB=zeros(length(SNRvec),1);
% 
% l=0;
% for L=SNRvec
%     l=l+1;
%     berFB(l)=sum(BER(l,n,idxA))/length(idxA);  
%     uppFB(l)=std(BER(l,n,idxA));
%     lowFB(l)=min(uppFB(l),berFB(l));
%     aderFB(l)=sum(pOnOff(l,n,idxA))/length(idxA);
% end

figure,%errorbar(SNRvec,ber,low,upp,'r');
plot(SNRvec,ber);
ylabel('BER'),xlabel('-10 log(\sigma_n^2)')
figure,plot(SNRvec,ader);
ylabel('ADER'),xlabel('-10 log(\sigma_n^2)')
figure,plot(SNRvec, dep);
legend('IFU-FSM'), ylabel('DEP'),xlabel('-10 log(\sigma_n^2)')
figure,plot(SNRvec, udep);
legend('IFU-FSM'), ylabel('UDEP'),xlabel('-10 log(\sigma_n^2)')
figure,plot(SNRvec, mse);
legend('IFU-FSM'), ylabel('MSE'),xlabel('-10 log(\sigma_n^2)')

figure,plot(SNRvec, mem1,SNRvec, mem2,SNRvec, mem3);
legend('Mem1','Mem2','Mem3'), ylabel('Memorias'),xlabel('-10 log(\sigma_n^2)')


