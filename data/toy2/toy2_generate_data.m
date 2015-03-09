clear all;
load ../toy1/Nt5_SNR50.mat

Nt = 5;
param.gen.s2n = 0;
param.gen.L_true = 2*ones(1,Nt);
param.gen.varH = 1;
param.gen.Nt = Nt;
param.gen.burstLength = 150*ones(1,Nt);
param.gen.burstLengthStdFactor = 5;
param.gen.symbol0 = 0;

param.L = 2;
param.flag0 = 1;
M = 4;
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.Nr = 25;
param.T = 500;

Lnew = max(param.gen.L_true)-1;
newH = sqrt(0.5*param.gen.varH)*randn(param.Nr,param.gen.Nt,Lnew)+1i*sqrt(0.5*param.gen.varH)*randn(param.Nr,param.gen.Nt,Lnew);
data.channel = cat(3,data.channel,newH);

dataIni = data;
clear data;

figure, plot(dataIni.seq')

for Nt=2:5
    for SNR=[-25:5:10 20 50]
        s2n = 10^(-SNR/10);
        
        data = dataIni;
        data.channel = data.channel(:,1:Nt,:);
        data.seq = data.seq(1:Nt,:);
        data.symbols = data.symbols(1:Nt,:);
        
        noise = sqrt(s2n/2)*randn(param.Nr,param.T)+1i*sqrt(s2n/2)*randn(param.Nr,param.T);
        
        data.obs = zeros(param.Nr,param.T);
        for ll=1:max(param.gen.L_true)
            data.obs = data.obs+data.channel(:,:,ll)*[zeros(Nt,ll-1) data.symbols(:,1:param.T-ll+1)];
        end
        data.obs = data.obs+noise;
        %figure; plot(data.seq');
        
        save(['Nt' num2str(Nt) '_SNR' num2str(SNR) '.mat'],'data');
    end
end
