close all;clear all;
Nvec=[2 3 5];
nvec=1:3;
load('resultadosMIMO_SNR.mat');


SNRvec=-5:5;
Memmean=zeros(3,length(SNRvec),5);
Nmean=zeros(3,length(SNRvec));
dep=zeros(3,length(SNRvec));
udep=zeros(3,length(SNRvec));

for n=nvec
    N=Nvec(n);
    l=0;
    for L=SNRvec
        l=l+1;
        idxA=find(Nfin(l,n,:)==N & Mfin(l,n,:)==1);
        dep(n,l)=1-sum(Mfin(l,n,:) & Nfin(l,n,:)==N)/100;
        udep(n,l)=1-sum(Nfin(l,n,:)==N)/100;
        Nmean(n,l)=mean(Nfin(l,n,:));
        
        idxA=find(Nfin(l,n,:)==N);
        Memmean(n,l,:)=mean(MEM(l,n,idxA,:),3);
    end
end

figure, h=plot(SNRvec,dep,'linewidth',2)
set(gca,'fontsize',14)
set(h(1),'marker','o', 'linestyle', '-.','markersize',8);
set(h(2),'marker','+', 'linestyle', '--','markersize',10);
set(h(3),'marker','x', 'linestyle', '-','markersize',10);
ylabel('DEP'),xlabel('-10 log(\sigma_n^2)')
grid on
legend('N_t=2','N_t=3','N_t=5');
set(gca,'xtick',SNRvec)
% figurapdf(6,4)
% print -dpdf DEP.pdf


figure, h=plot(SNRvec,udep,'linewidth',2)
set(gca,'fontsize',14)
set(h(1),'marker','o', 'linestyle', '-.','markersize',8);
set(h(2),'marker','+', 'linestyle', '--','markersize',10);
set(h(3),'marker','x', 'linestyle', '-','markersize',10);

ylabel('UDEP'),xlabel('-10 log(\sigma_n^2)')
grid on
legend('N_t=2','N_t=3','N_t=5');
set(gca,'xtick',SNRvec)
figurapdf(6,4)
print -dpdf UDEP.pdf


figure, h=plot(SNRvec,Nmean,'linewidth',2)
set(gca,'fontsize',14)
set(h(1),'marker','o', 'linestyle', '-.','markersize',8);
set(h(2),'marker','+', 'linestyle', '--','markersize',10);
set(h(3),'marker','x', 'linestyle', '-','markersize',10);
ylabel('# Transmitters'),xlabel('-10 log(\sigma_n^2)')
grid on
legend('N_t=2','N_t=3','N_t=5','location','northwest');
set(gca,'xtick',SNRvec)

figurapdf(6,4)
print -dpdf Nmean.pdf


Txs=['Tx ';
    'Tx ';
    'Tx ';
    'Tx ';
    'Tx ';
    'Tx '];
rs=[', r_1^{}=';
    ', r_2^{}=';
    ', r_3^{}=';
    ', r_4^{}=';
    ', r_5^{}=';
    ', r_6^{}='];

markers=['x','+','o','^','s'];
for n=nvec
    if n==1
        Mvec=[2 2]-1;
    elseif n==2
        Mvec=[3 2 1]-1;
    else 
        Mvec=[4 4 3 2 2]-1;
    end
    figure,h=plot(SNRvec,reshape(Memmean(n,:,1:Nvec(n)),length(SNRvec),Nvec(n)),'linewidth',2);
    set(gca,'fontsize',14)
    for i=1:Nvec(n)
        if i<3
            set(h(i),'marker',markers(i), 'linestyle', '-','markersize',10);
        else
            set(h(i),'marker',markers(i), 'linestyle', '-','markersize',8);
        end
    end
    ylabel('Memory error'),xlabel('-10 log(\sigma_n^2)')
    set(gca,'xtick',SNRvec)
    grid on
    legend([Txs(1:Nvec(n),:) num2str([1:Nvec(n)]') rs(1:Nvec(n),:) num2str(Mvec')],'location','southwest');
    figurapdf(6,4)
    print ('-dpdf', ['Mem' num2str(n) '.pdf']);
    

end