%% Simulations for 2 transmitters and channel lengh 1
Lvec=[10];%Number of receivers
SNRvec=-3:1:5; %SNR

%% Parameteres to be set
% Number of receiver indixes
L1=1;



%% Plotting
load ('resultPruebas2.mat')
%close all
f=0;
fChi=0;
% %Q
% figure,plot(-10:10,sum(Qfin{Qini}(11,:,:)==3,3)/100, -10:10,sum(Qfin{Qini}(6,:,:)==3,3)/100, -10:10,sum(Qfin{Qini}(2,:,:)==3,3)/100)
% 
% %M
% plot(-10:10,sum(Mfin{Qini}(11,:,:)==3,3)/100, -10:10,sum(Mfin{Qini}(6,:,:)==3,3)/100, -10:10,sum(Mfin{Qini}(2,:,:)==3,3)/100)


%T1_DEP.eps
figure,plot(SNRvec,1-reshape(sum(acierto(f+1,fChi+1,L1,:,:),5),1,length(SNRvec))/100,SNRvec,1-reshape(sum(aciertoM(f+1,fChi+1,L1,:,:),5),1,length(SNRvec))/100,'linewidth',2)
grid on
xlabel('SNR (dB)'); ylabel('DEP')
legend(['N_r=' num2str(Lvec(L1))],'Location','SouthWest');
% axis([-8 8 0 1]);
% figurapdf(5.5,2.5);
% print -dpdf T1_DEP.pdf

