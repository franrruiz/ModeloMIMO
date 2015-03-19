function simClusterSyn(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,itCluster,simId)

addpath(genpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab'));

randn('seed',round(sum(1e5*clock)+itCluster));
rand('seed',round(sum(1e5*clock)+itCluster));

saveFolder = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_Ltrue' num2str(Ltrue) '_L' num2str(L) '_SNR' num2str(SNR) '_lHead' num2str(lHead), '_onOff' num2str(onOffModel)];
saveFile = [saveFolder '/itCluster' num2str(itCluster)];

if(~isdir(saveFolder))
    mkdir(saveFolder);
end
if(~isdir(saveFile))
    mkdir(saveFile);
end

if(exist([saveFile '.mat'],'file'))
    return;
end

%% Configuration parameters
param.Nr = Nr;                        % Number of receivers
param.T  = T;                         % Length of the sequence
% M for the M-QAM constellation
M = 2^M;
param.constellation = qammod(0:M-1,M,[],'gray');
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = L;        % Channel memory to be considered during inference
param.Niter = Niter;  % Number of iterations of the sampler
param.saveCycle = 200;
param.storeIters = 2000;
param.header = ones(1,lHead);
param.onOffModel = onOffModel;

%% Generate data
if(simId>=3)
    noiseVar = 10^(-SNR/10);
else
    noiseVar = param.Nr*10^(-SNR/10);
end
if(log2(M)==1)
    noiseVar = 2*noiseVar;
end
param.gen.s2n = noiseVar;
param.gen.L_true = Ltrue*ones(1,Nt);
param.gen.varH = 1;
if(log2(M)==1)
    param.gen.varH = 2*param.gen.varH;
end
param.gen.Nt = Nt;
param.gen.burstLength = round(T/2)*ones(1,param.gen.Nt);
param.gen.burstLengthStdFactor = inf;
param.gen.symbol0 = 0;
param.gen.sparsityH = 0;

if(simId<=10)
    data = generate_data_bursts(param);
elseif(simId==11)
    % If simId==11, load data from file
    load(['/export/clusterdata/franrruiz87/ModeloMIMO/data/syn/' num2str(simId) '/T' num2str(param.T) '_itCluster' num2str(itCluster) '.mat']);
    data.channel = channelGen(1:param.Nr,1:param.gen.Nt,1:max(param.gen.L_true));
    data.symbols = symbolsGen{log2(M)}(1:param.gen.Nt,1:param.T);
    data.seq = seqGen{log2(M)}(1:param.gen.Nt,1:param.T);
    data.obs = zeros(param.Nr,param.T);
    
    noise = sqrt(param.gen.s2n/2)*noiseGen(1:param.Nr,1:param.T);
    if(log2(M)==1)
        data.symbols = real(data.symbols);
        data.channel = sqrt(2)*real(data.channel);
        noise = sqrt(2)*real(noise);
    end
    for t=1:param.T
        for i=0:max(param.gen.L_true)-1
            if(t-i<=0)
                data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*param.gen.symbol0*ones(param.gen.Nt,1);
            else
                data.obs(:,t) = data.obs(:,t) + data.channel(:,:,i+1)*data.symbols(:,t-i);
            end
        end
    end
    data.obs = data.obs+noise;
else
    error('I do not know how to generate data');
end

%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;
param.pgas.maxM = 40;
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
param.ep.eps = 5e-7;
param.ep.beta = 0.2;
param.ep.Niter = 15;
param.colGibbs.Niter = 1;
param.ffbs.Niter = 1;

%% Configuration parameters for BNP and inference method
param.infer.symbolMethod = 'pgas';
param.infer.sampleNoiseVar = 0;
param.infer.sampleChannel = 1;
param.infer.sampleVarH = 1;
param.bnp.betaSlice1 = 0.5;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;
param.bnp.Mini = 1;

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

%% Check if there are temporary files to be loaded
flagRecovered = 0;
itInit = 0;
it = param.saveCycle;
while(it<=param.Niter)
    if(exist([saveFile '/it' num2str(it) '.mat'],'file'))
        try
            % Try to load the temporary file
            load([saveFile '/it' num2str(it) '.mat']);
            % If success, then save current iteration and activate flag
            itInit = it;
            flagRecovered = 1;
            % Delete previous file (it-saveCycle) in order not to exceed disk quota
            if(exist([saveFile '/it' num2str(it-param.saveCycle) '.mat'],'file'))
                delete([saveFile '/it' num2str(it-param.saveCycle) '.mat']);
            end
        catch e
            % If the file exists but it is corrupt
            % (it happens sometimes when using the cluster machines)
            delete([saveFile '/it' num2str(it) '.mat']);
            % Set it so that next iteration is it-param.saveCycle
            it = it-2*param.saveCycle;
            itInit = 0;
            flagRecovered = 0;
        end
    end
    it = it+param.saveCycle;
end

%% Initialization
if(~flagRecovered)
    init.H = zeros(param.Nr,param.bnp.Mini,param.L);
    for ll=1:param.L
        init.H(:,:,ll) = sqrt(hyper.s2h*exp(-hyper.lambda*(ll-1)))*(randn(param.Nr,param.bnp.Mini,1)+1i*randn(param.Nr,param.bnp.Mini,1));
    end
    if(~param.infer.sampleNoiseVar)
        init.s2y = noiseVar;      % INITIALIZE s2y TO THE GROUND TRUTH
    else
        init.s2y = 20*rand(1);    % INITIALIZE s2y TO SOME LARGE VALUE
    end
    init.s2H = hyper.s2h*exp(-hyper.lambda*(0:param.L-1));  % INITIALIZE s2H TO ITS MEAN VALUE
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
    ADER = zeros(1,param.Niter+1);
    SER_ALL = zeros(1,param.Niter+1);
    SER_ACT = zeros(1,param.Niter+1);
    MMSE = zeros(1,param.Niter+1);
    LLH = zeros(1,param.Niter+1);
    M_EST = zeros(1,param.Niter+1);
    samplesAll = cell(1,param.storeIters);
end
for it=itInit+1:param.Niter
    %% Algorithm
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % For PGAS, check that the number of current chains does not exceed maxM
    if(strcmp(param.infer.symbolMethod,'pgas'))
        if(size(samples.seq,1)>param.pgas.maxM)
            param.pgas.maxM = size(samples.seq,1);
            param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
        end
    end
    
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
    samples.s2y = sample_post_s2y(data,samples,hyper,param);
    % -Sample the variance of the channel coefficients
    samples.s2H = sample_post_s2H(data,samples,hyper,param);
    
    %% Store current sample
    if(it>param.Niter-param.storeIters)
        samplesAll{it-param.Niter+param.storeIters} = samples;
    end
    
    %% Evaluation
    % Trace of the estimated number of transmitters
    M_EST(it) = sum(sum(samples.seq~=0,2)>0);
    % Trace of the log-likelihood
    LLH(it) = compute_llh(data,samples,hyper,param);
    % SER, ADER
    if(it==param.Niter)
        [ADER(it) SER_ALL(it) SER_ACT(it) MMSE(it) vec_ord rot] = compute_error_rates(data,samples,hyper,param);
    end
    
    %% Save temporary result file
    if(mod(it,param.saveCycle)==0)
        save([saveFile '/it' num2str(it) '.mat'],'data','init','samples','ADER','SER_ALL','SER_ACT','MMSE','LLH','M_EST','samplesAll');
        % If successfully saved, delete previous temporary file
        if(exist([saveFile '/it' num2str(it-param.saveCycle) '.mat'],'file'))
            delete([saveFile '/it' num2str(it-param.saveCycle) '.mat']);
        end
    end
end

%% Final evaluation of performance
Zaux = zeros(size(samples.Z,1),param.T,1+length(param.constellation));
auxSample.s2H = zeros(1,param.L);
for it=1:param.storeIters
    if(size(samplesAll{it}.seq,1)>size(Zaux,1))
        Zaux = cat(1,Zaux,zeros(size(samplesAll{it}.seq,1)-size(Zaux,1),param.T,1+length(param.constellation)));
    end
    for t=1:param.T
        for m=1:size(samplesAll{it}.seq,1)
            Zaux(m,t,samplesAll{it}.seq(m,t)+1) = 1+Zaux(m,t,samplesAll{it}.seq(m,t)+1);
        end
    end
    auxSample.s2H = auxSample.s2H+(samplesAll{it}.s2H/param.storeIters);
end
[valnul auxIdx] = max(Zaux,[],3);
auxConstellation = [0 param.constellation];
auxSample.seq = auxIdx-1;
auxSample.Z = auxConstellation(auxIdx);
auxSample.s2y = samples.s2y;
[valnul auxSample.H] = sample_post_H(data,auxSample,hyper,param);

[ADER(param.Niter+1) SER_ALL(param.Niter+1) SER_ACT(param.Niter+1) MMSE(param.Niter+1) vec_ord rot ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates(data,auxSample,hyper,param);
LLH(param.Niter+1) = compute_llh(data,auxSample,hyper,param);
M_EST(param.Niter+1) = sum(sum(auxSample.seq~=0,2)>0);

%% Save result
save([saveFile '.mat'],'data','init','samples','ADER','SER_ALL','SER_ACT','MMSE','LLH','M_EST','samplesAll','*_indiv');

% If successfully saved, detele previous temporary file
if(exist([saveFile '/it' num2str(param.Niter) '.mat'],'file'))
    delete([saveFile '/it' num2str(param.Niter) '.mat']);
end

