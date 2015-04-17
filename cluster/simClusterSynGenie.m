function simClusterSynGenie(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,Nparticles,flagParallel,itCluster,simId)
% 
% Initially, it created variables ADER_PGAS, ADER_BCJR, etc.
% It was computed with am=0.95, bm=0.05
% Now, it creates variables ADER_PGAS2, ADER_BCJR2, etc.
% It is computed with am=0.995, bm=0.005
% 
% Furthermore, initially it returned 4000 samples for the PGAS (40%), now it
% returns only 2000 (20%)
% 
% More importantly, variables ending in "2" had been computed in a wrong
% way, they should be ignored and replaced with variables ending in "3"
% 
% Variables ending in "4" have am=499/500 and bm=1/500
% 

addpath(genpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab'));

randn('seed',round(sum(1e5*clock)+itCluster));
rand('seed',round(sum(1e5*clock)+itCluster));

saveFolder = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId)...
              '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_Ltrue' num2str(Ltrue) '_L' num2str(L) '_SNR' num2str(SNR) '_lHead' num2str(lHead), '_onOff' num2str(onOffModel) '_Npart' num2str(Nparticles)];
saveFile = [saveFolder '/itCluster' num2str(itCluster)];
saveTmpFolder = [saveFile '/pgas'];

% Exit if BNP+PGAS has not finished yet
if(~exist([saveFile '.mat'],'file'))
    return;
end
% Load data
load([saveFile '.mat'],'data','ADER_PGAS4');
% Exit if this simulation already run
if(exist('ADER_PGAS4','var'))
    return;
end
if(~isdir(saveTmpFolder))
    mkdir(saveTmpFolder);
end

%% Configuration parameters
param.Nr = Nr;                        % Number of receivers
param.T  = T;                         % Length of the sequence
% M for the M-QAM constellation
M = 2^M;
param.constellation = qammod(0:M-1,M,[],'gray');
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = Ltrue;        % Channel memory to be considered during inference
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = 200;
param.storeIters = round(Niter/5);
param.header = [];
param.onOffModel = 0;

%% Choose noiseVar
noiseVar = 10^(-SNR/10);
if(log2(M)==1)
    noiseVar = 2*noiseVar;
end

%% Check if there are temporary files to be loaded
flagRecovered = 0;
itInit = 0;
it = param.saveCycle;
while(it<=param.Niter)
    if(exist([saveTmpFolder '/it' num2str(it) '.mat'],'file'))
        try
            % Try to load the temporary file
            load([saveTmpFolder '/it' num2str(it) '.mat']);
            % If success, then save current iteration and activate flag
            itInit = it;
            flagRecovered = 1;
            % Delete previous file (it-saveCycle) in order not to exceed disk quota
            if(exist([saveTmpFolder '/it' num2str(it-param.saveCycle) '.mat'],'file'))
                delete([saveTmpFolder '/it' num2str(it-param.saveCycle) '.mat']);
            end
        catch e
            % If the file exists but it is corrupt
            % (it happens sometimes when using the cluster machines)
            delete([saveTmpFolder '/it' num2str(it) '.mat']);
            % Set it so that next iteration is it-param.saveCycle
            it = it-2*param.saveCycle;
            itInit = 0;
            flagRecovered = 0;
        end
    end
    it = it+param.saveCycle;
end

%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = (0.5*param.T-1)/(0.5*param.T);
param.bcjr.p2 = 1/(0.5*param.T);
param.pgas.N_PF = Nparticles;
param.pgas.N_PG = Nparticles;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = size(data.symbols,1);
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.pgas.flagParallel = flagParallel;  %%% IF 1, IT CAN BE RUN ONLY OVER MULTICORE MACHINES
param.ep.eps = 5e-7;
param.ep.beta = 0.2;
param.ep.Niter = 15;
param.colGibbs.Niter = 1;
param.ffbs.Niter = 1;

%% Configuration parameters for BNP and inference method
param.infer.sampleNoiseVar = 0;
param.infer.sampleChannel = 0;
param.infer.sampleVarH = 0;
param.infer.simulatedTempering = 0;
param.infer.addArtificialNoise = 0;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = size(data.symbols,1);

%% Hyperparameters
hyper.s2h = 1;      % E[s2H(r)]=s2h*exp(-lambda*(r-1))
if(M==1)
    hyper.s2h = 2*hyper.s2h;
end
hyper.lambda = 0.5; % E[s2H(r)]=s2h*exp(-lambda*(r-1))
hyper.kappa = 1;    % Std[s2H(r)]=kappa*E[s2H(r)]
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization
samplesPGAS.H = data.channel;    % INITIALIZE channel TO THE GROUND TRUTH
samplesPGAS.s2y = noiseVar;      % INITIALIZE s2y TO THE GROUND TRUTH
samplesPGAS.s2H = hyper.s2h*exp(-hyper.lambda*(0:param.L-1));  % INITIALIZE s2H TO ITS MEAN VALUE
samplesPGAS.am = ((0.5*param.T-1)/(0.5*param.T))*ones(param.bnp.Mini,1);
samplesPGAS.bm = (1/(0.5*param.T))*ones(param.bnp.Mini,1);
samplesPGAS.Z = zeros(param.bnp.Mini,param.T);
samplesPGAS.seq = zeros(param.bnp.Mini,param.T);
samplesPGAS.nest = zeros(2,2,param.bnp.Mini);
samplesPGAS.nest(1,1,:) = param.T;
samplesPGAS.slice = 0;
samplesPGAS.epAcc = 0;

samplesFFBS = samplesPGAS;   % Initialize FFBS as PGAS

%% Inference using PGAS and FFBS
if(~flagRecovered)
    auxSamplePGAS = samplesPGAS;
    ZauxPGAS = zeros(param.bnp.Mini,param.T,1+length(param.constellation));
    auxSampleFFBS = samplesFFBS;
    ZauxFFBS = zeros(param.bnp.Mini,param.T,1+length(param.constellation));
end
for it=itInit+1:param.Niter
    %% Run PGAS
    [auxSamplePGAS.Z auxSamplePGAS.seq] = pgas_main(data,auxSamplePGAS,hyper,param);
    % Update ZauxPGAS (to average)
    if(it>param.Niter-param.storeIters)
        for t=1:param.T
            for m=1:size(auxSamplePGAS.Z,1)
                ZauxPGAS(m,t,auxSamplePGAS.seq(m,t)+1) = 1+ZauxPGAS(m,t,auxSamplePGAS.seq(m,t)+1);
            end
        end
    end
    %% Run FFBS
    if((length(param.constellation)+1)^(2*Ltrue)<1e6)
        [auxSampleFFBS.Z auxSampleFFBS.seq valnul] = sample_Z_FFBS(data,auxSampleFFBS,hyper,param);
        % Update ZauxFFBS (to average)
        if(it>param.Niter-param.storeIters)
            for t=1:param.T
                for m=1:size(auxSampleFFBS.Z,1)
                    ZauxFFBS(m,t,auxSampleFFBS.seq(m,t)+1) = 1+ZauxFFBS(m,t,auxSampleFFBS.seq(m,t)+1);
                end
            end
        end
    end
    %% Save tmp results to file
    if(mod(it,param.saveCycle)==0)
        save([saveTmpFolder '/it' num2str(it) '.mat'],'auxSamplePGAS','ZauxPGAS','auxSampleFFBS','ZauxFFBS');
        % If successfully saved, delete previous temporary file
        if(exist([saveTmpFolder '/it' num2str(it-param.saveCycle) '.mat'],'file'))
            delete([saveTmpFolder '/it' num2str(it-param.saveCycle) '.mat']);
        end
    end
end

%% Performance of PGAS
[valnul auxIdx] = max(ZauxPGAS,[],3);
auxConstellation = [0 param.constellation];
auxSamplePGAS.seq = auxIdx-1;
auxSamplePGAS.Z = auxConstellation(auxIdx);
[ADER_PGAS4 SER_ALL_PGAS4 SER_ACT_PGAS4 MMSE_PGAS4 vec_ord rot ADER_PGAS4_indiv SER_ALL_PGAS4_indiv SER_ACT_PGAS4_indiv MMSE_PGAS4_indiv] = ...
    compute_error_rates_genie(data,auxSamplePGAS,hyper,param);

%% Performance of FFBS
ADER_FFBS4 = NaN;
SER_ALL_FFBS4 = NaN;
SER_ACT_FFBS4 = NaN;
MMSE_FFBS4 = NaN;
ADER_FFBS4_indiv = NaN;
SER_ALL_FFBS4_indiv = NaN;
SER_ACT_FFBS4_indiv = NaN;
MMSE_FFBS4_indiv = NaN;
if((length(param.constellation)+1)^(2*Ltrue)<1e6)
    [valnul auxIdx] = max(ZauxFFBS,[],3);
    auxConstellation = [0 param.constellation];
    auxSampleFFBS.seq = auxIdx-1;
    auxSampleFFBS.Z = auxConstellation(auxIdx);
    [ADER_FFBS4 SER_ALL_FFBS4 SER_ACT_FFBS4 MMSE_FFBS4 vec_ord rot ADER_FFBS4_indiv SER_ALL_FFBS4_indiv SER_ACT_FFBS4_indiv MMSE_FFBS4_indiv] = ...
        compute_error_rates_genie(data,auxSampleFFBS,hyper,param);
end

%% Inference using BCJR
ADER_BCJR4 = NaN;
SER_ALL_BCJR4 = NaN;
SER_ACT_BCJR4 = NaN;
MMSE_BCJR4 = NaN;
ADER_BCJR4_indiv = NaN;
SER_ALL_BCJR4_indiv = NaN;
SER_ACT_BCJR4_indiv = NaN;
MMSE_BCJR4_indiv = NaN;
if((length(auxConstellation)^(2*param.L*param.bnp.Mini)<1e6))
    auxSample = samplesPGAS;
    [auxSample.Z qt_red Simb_red] = bcjr_main(data,samplesPGAS,hyper,param);
    auxSample.seq = zeros(size(auxSample.Z));
    for q=1:length(auxConstellation)
        idx = (abs(auxSample.Z-auxConstellation(q))<min(abs(param.constellation))/10);
        auxSample.seq(idx) = q-1;
    end
    [ADER_BCJR4 SER_ALL_BCJR4 SER_ACT_BCJR4 MMSE_BCJR4 vec_ord rot ADER_BCJR4_indiv SER_ALL_BCJR4_indiv SER_ACT_BCJR4_indiv MMSE_BCJR4_indiv] = ...
        compute_error_rates_genie(data,auxSample,hyper,param);
end

%% Save results
save([saveFile '.mat'],'ADER_PGAS4','SER_ALL_PGAS4','SER_ACT_PGAS4','MMSE_PGAS4',...
                       'ADER_FFBS4','SER_ALL_FFBS4','SER_ACT_FFBS4','MMSE_FFBS4',...
                       'ADER_BCJR4','SER_ALL_BCJR4','SER_ACT_BCJR4','MMSE_BCJR4',...
                       '*4_indiv',...
                       '-append');
                   
%% If successfully saved, detele previous temporary file
if(exist([saveTmpFolder '/it' num2str(param.Niter) '.mat'],'file'))
    delete([saveTmpFolder '/it' num2str(param.Niter) '.mat']);
end
