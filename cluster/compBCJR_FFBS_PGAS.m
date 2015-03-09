function compBCJR_FFBS_PGAS(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,itCluster,simId)

addpath(genpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab'));

randn('seed',round(sum(1e5*clock)+itCluster));
rand('seed',round(sum(1e5*clock)+itCluster));

saveFolder = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_Ltrue' num2str(Ltrue) '_L' num2str(L) '_SNR' num2str(SNR) '_lHead' num2str(lHead), '_onOff' num2str(onOffModel)];
saveFile = [saveFolder '/itCluster' num2str(itCluster)];

if(~isdir(saveFolder))
    mkdir(saveFolder);
end

if(exist([saveFile '.mat'],'file'))
    return;
end

%% Configuration parameters
noiseVar = 10^(-SNR/10);
if(M==1)
    noiseVar = 2*noiseVar;
end
param.Nr = Nr;                        % Number of receivers
param.T  = T;                         % Length of the sequence
% M for the M-QAM constellation
M = 2^M;
param.constellation = qammod(0:M-1,M,[],'gray');
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = L;        % Channel memory to be considered during inference
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = inf;
param.storeIters = round(Niter/2);
param.header = ones(1,lHead);
param.onOffModel = onOffModel;

%% Generate data
param.gen.Nt = Nt;   % Number of transmitters (to generate data)
param.gen.L_true = Ltrue*ones(1,Nt);  % Channel memory for each transmitter (including tap 0)
param.gen.symbol0 = 0;   % Symbols transmitted before (and after) the transmission
param.gen.s2n = noiseVar;                    % Variance of the noise (to generate data)
param.gen.varH = 1;                   % Variance of the channel (to generate data)
param.gen.burstLength = round(T/2)*ones(1,param.gen.Nt);  % Mean length of bursts (to generate data)
param.gen.burstLengthStdFactor = inf;
param.gen.sparsityH = 0;

% Generate data:
data = generate_data_bursts(param);

%% Configuration parameters for BCJR, FFBS and PGAS
param.bcjr.p1 = 0.995;
param.bcjr.p2 = 0.005;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = param.gen.Nt;
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.ffbs.Niter = 1;

%% Create struct samples
samples.H = data.channel;
samples.s2y = param.gen.s2n;
samples.am = param.bcjr.p1*ones(param.gen.Nt,1);
samples.bm = param.bcjr.p2*ones(param.gen.Nt,1);
samples.s2H = param.gen.varH*ones(1,param.L);
if(param.flag0)
    samples.seq = randint(param.gen.Nt,param.T,[0 M]);
    auxConstellation = [0 param.constellation];
    samples.Z = auxConstellation(samples.seq+1);
end

%% Initialization and inference
% (1) PGAS
samplesPGAS = samples;
P_PGAS = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
for it=1:param.Niter
    [samplesPGAS.Z samplesPGAS.seq] = pgas_main(data,samplesPGAS,[],param);
    if(it>param.storeIters)
        for t=1:param.T
            for n=1:param.gen.Nt
                P_PGAS(n,t,1) = P_PGAS(n,t,1)+sum(samplesPGAS.seq(n,t)==0)/param.storeIters;
                for q=1:length(param.constellation)
                    P_PGAS(n,t,q+1) = P_PGAS(n,t,q+1)+sum(samplesPGAS.seq(n,t)==q)/param.storeIters;
                end
            end
        end
    end
end

% (2) FFBS
samplesFFBS = samples;
P_FFBS = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
for it=1:param.Niter
    [samplesFFBS.Z samplesFFBS.seq valnul] = sample_Z_FFBS(data,samplesFFBS,[],param);
    if(it>param.storeIters)
        for t=1:param.T
            for n=1:param.gen.Nt
                P_FFBS(n,t,1) = P_FFBS(n,t,1)+sum(samplesFFBS.seq(n,t)==0)/param.storeIters;
                for q=1:length(param.constellation)
                    P_FFBS(n,t,q+1) = P_FFBS(n,t,q+1)+sum(samplesFFBS.seq(n,t)==q)/param.storeIters;
                end
            end
        end
    end
end

% (3) BCJR
samplesBCJR = samples;
P_BCJR = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
[Sest qt_red Simb_red] = bcjr_main(data,samplesBCJR,[],param);
for t=1:param.T
    for n=1:param.gen.Nt
        idx = find(Simb_red(:,n)==0);
        P_BCJR(n,t,1) = sum(qt_red(idx,t));
        for q=1:length(param.constellation)
            idx = find(Simb_red(:,n)==param.constellation(q));
            P_BCJR(n,t,q+1) = sum(qt_red(idx,t));
        end
    end
end

%% Save results
save([saveFile '.mat'],'P_PGAS','P_FFBS','P_BCJR','samplesPGAS','samplesFFBS','data','samples');

