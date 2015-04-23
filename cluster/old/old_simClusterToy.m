function simClusterToy(toyType,subToyType,T,Nt,Nr,M,L,SNR,Niter,itCluster)

addpath(genpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab'));

randn('seed',round(sum(1e5*clock)+itCluster));
rand('seed',round(sum(1e5*clock)+itCluster));

if(toyType>=1)
    loadFile = ['/export/clusterdata/franrruiz87/ModeloMIMO/data/toy' num2str(toyType) '/Nt' num2str(Nt) '_SNR' num2str(SNR) '.mat'];
end

saveFolder = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/toy' num2str(toyType) '/sub' num2str(subToyType)];
saveFile = [saveFolder '/Nt' num2str(Nt) '_SNR' num2str(SNR) '_itCluster' num2str(itCluster)];
if(subToyType==1)
    param.infer.sampleNoiseVar = 0;  % subTypes 1
else
    param.infer.sampleNoiseVar = 1;  % subTypes 2-4
end
param.infer.symbolMethod = 'pgas';

if(~isdir(saveFolder))
    mkdir(saveFolder);
end
if(~isdir(saveFile))
    mkdir(saveFile);
end

if(exist([saveFile '.mat'],'file'))
    return;
end

%% Load data
load(loadFile,'data');

%% Configuration parameters
param.Nr = Nr;                        % Number of receivers
param.T  = T;                         % Length of the sequence
% M for the M-QAM constellation
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = L;        % Channel memory to be considered during inference
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = 1000;

%% Configuration parameters for BCJR, PGAS, EP and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.ep.eps = 5e-7;
param.ep.beta = 0.2;
param.ep.Niter = 15;
param.colGibbs.Niter = 1;

%% Configuration parameters for BNP
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;

%% Hyperparameters
hyper.s2h = 1;      % Variance of H_r ~ CN(0,s2h*exp(-lambda*r)*I)
hyper.lambda = 0.2; % Variance of H_r ~ CN(0,s2h*exp(-lambda*r)*I)
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 0.1; % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 2;   % Parameter for bm ~ Beta(gamma1,gamma2)
if(subToyType==3)
    hyper.tau = 11/9;      % Parameter for s2y ~ IG(tau,nu)
    hyper.nu = 20/9;       % Parameter for s2y ~ IG(tau,nu)
else
    hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
    hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)
end

%% Check if there are temporary files to be loaded
flagRecovered = 0;
itInit = 0;
for it=param.saveCycle:param.saveCycle:param.Niter
    if(exist([saveFile '/it' num2str(it) '.mat'],'file'))
        load([saveFile '/it' num2str(it) '.mat']);
        itInit = it;
        flagRecovered = 1;
    end
end

%% Initialization
if(~flagRecovered)
    init.H = zeros(param.Nr,param.bnp.Mini,param.L);
    for ll=1:param.L
        init.H(:,:,ll) = sqrt(hyper.s2h*exp(-hyper.lambda*(ll-1)))*(randn(param.Nr,param.bnp.Mini,1)+1i*randn(param.Nr,param.bnp.Mini,1));
    end
    if(~param.infer.sampleNoiseVar)
        init.s2y = 10^(-SNR/10);  % INITIALIZE s2y TO THE GROUND TRUTH
    else
        init.s2y = 20*rand(1);    % INITIALIZE s2y TO SOME LARGE VALUE
    end
    init.am = 0.95*ones(param.bnp.Mini,1);
    init.bm = 0.05*ones(param.bnp.Mini,1);
    init.Z = zeros(param.bnp.Mini,param.T);
    init.seq = zeros(param.bnp.Mini,param.T);
    init.nest = zeros(2,2,param.bnp.Mini);
    init.nest(1,1,:) = param.T;
    init.slice = 0;
    init.epAcc = 0;
    samples = init;
end

%% Inference
if(~flagRecovered)
    ADER = zeros(1,param.Niter);
    SER_ALL = zeros(1,param.Niter);
    SER_ACT = zeros(1,param.Niter);
    LLH = zeros(1,param.Niter);
    M_EST = zeros(1,param.Niter);
end
for it=itInit+1:param.Niter
    %% Algorithm
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,hyper,param);
    % -Compute some statistics of interest
    if(strcmp(param.infer.symbolMethod,'pgas'))
        
    elseif(strcmp(param.infer.symbolMethod,'ep'))
        samples.epAcc = samples.epAcc+out;
    end
    
    % Step 3)
    % -Remove unused chains
    samples = sample_remove_unused(data,samples,hyper,param);
    
    % Step 4)
    % -Sample the transition probabilities (semi-ordered construction)
    [samples.am samples.bm]= sample_post_transitionProb(data,samples,hyper,param);
    
    % Step 5)
    % -Sample the channel H
    samples.H = sample_post_H(data,samples,hyper,param);
    % -Sample the noise variance
    if(subToyType~=4 || it>100)
        samples.s2y = sample_post_s2y(data,samples,hyper,param);
    end
    %% Evaluation
    % Trace of the estimated number of transmitters
    M_EST(it) = sum(sum(samples.seq~=0,2)>0);
    % Trace of the log-likelihood
    LLH(it) = compute_llh(data,samples,hyper,param);
    % SER, ADER
    if((mod(it,200)==0)&&(size(samples.Z,1)<=11))
        [ADER(it) SER_ALL(it) SER_ACT(it) vec_ord rot] = compute_error_rates(data,samples,hyper,param);
    end
    
    if(mod(it,param.saveCycle)==0)
        save([saveFile '/it' num2str(it) '.mat'],'init','samples','ADER','SER_ALL','SER_ACT','LLH','M_EST');
    end
end

save([saveFile '.mat'],'init','samples','ADER','SER_ALL','SER_ACT','LLH','M_EST');

