
Nt = 5;
param.gen.s2n = 0;
param.gen.L_true = ones(1,Nt);
param.gen.varH = 1;
param.gen.Nt = Nt;
param.gen.burstLength = 150*ones(1,Nt);
param.gen.burstLengthStdFactor = 5;
param.gen.symbol0 = 0;

param.L = 1;
param.flag0 = 1;
M = 4;
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.Nr = 25;
param.T = 500;

dataIni = generate_data_bursts(param);
figure, plot(dataIni.seq')

idxOrder = randperm(5);
for Nt=2:5
    for SNR=[-10 0 10 20 50]
        s2n = 10^(-SNR/10);
        
        data = dataIni;
        data.channel = data.channel(:,idxOrder(1:Nt));
        data.seq = data.seq(idxOrder(1:Nt),:);
        data.symbols = data.symbols(idxOrder(1:Nt),:);
        noise = sqrt(s2n/2)*randn(param.Nr,param.T)+1i*sqrt(s2n/2)*randn(param.Nr,param.T);
        
        data.obs = data.channel*data.symbols+noise;

        figure; plot(data.seq');
        
        save(['Nt' num2str(Nt) '_SNR' num2str(SNR) '.mat'],'data');
    end
end



