function simClusterSynAll(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,itCluster,simId)

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
if(M==1)
    noiseVar = 2*noiseVar;
end
param.gen.s2n = noiseVar;
param.gen.L_true = Ltrue*ones(1,Nt);
param.gen.varH = 1;
if(M==1)
    param.gen.varH = 2*param.gen.varH;
end
param.gen.Nt = Nt;
param.gen.burstLength = round(T/2)*ones(1,param.gen.Nt);
param.gen.burstLengthStdFactor = inf;
param.gen.symbol0 = 0;
param.gen.sparsityH = 0.1;

data = generate_data_bursts(param);

%% Configuration parameters for BCJR, PGAS, EP, FFBS and collapsed Gibbs
param.bcjr.p1 = 0.995;
param.bcjr.p2 = 0.005;
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

%% Algorithms to be run
flag_genieFFBS = 0;
if((length(param.constellation)+1)^(2*Ltrue)<1e6)
    flag_genieFFBS = 1;
end
flag_genieBCJR = 0;
if(((1+length(param.constellation))^(2*Ltrue*Nt)<1e6)&&(Nt<=4))
    flag_genieBCJR = 1;
end

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
            % Set 'it' so that next iteration is it-param.saveCycle
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
    
    initGenie.H = data.channel;
    initGenie.s2y = noiseVar;
    initGenie.s2H = param.gen.varH*ones(1,Ltrue);
    initGenie.am = param.bcjr.p1*ones(Nt,1);
    initGenie.bm = param.bcjr.p2*ones(Nt,1);
    initGenie.Z = zeros(Nt,param.T);
    initGenie.seq = zeros(Nt,param.T);
    initGenie.nest = zeros(2,2,Nt);
    initGenie.nest(1,1,:) = param.T;
    initGenie.slice = 0;
    initGenie.epAcc = 0;
    
    samples_genieFFBS = initGenie;
    samples_geniePGAS = initGenie;
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

    ADER_genieFFBS = zeros(1,param.Niter+1);
    SER_ALL_genieFFBS = zeros(1,param.Niter+1);
    SER_ACT_genieFFBS = zeros(1,param.Niter+1);
    LLH_genieFFBS = zeros(1,param.Niter+1);
    samplesAll_genieFFBS = cell(1,param.storeIters);

    ADER_geniePGAS = zeros(1,param.Niter+1);
    SER_ALL_geniePGAS = zeros(1,param.Niter+1);
    SER_ACT_geniePGAS = zeros(1,param.Niter+1);
    LLH_geniePGAS = zeros(1,param.Niter+1);
    samplesAll_geniePGAS = cell(1,param.storeIters);

    ADER_genieBCJR = zeros(1,param.Niter+1);
    SER_ALL_genieBCJR = zeros(1,param.Niter+1);
    SER_ACT_genieBCJR = zeros(1,param.Niter+1);
    LLH_genieBCJR = zeros(1,param.Niter+1);
end

for it=itInit+1:param.Niter
    %% BNP (PGAS Algorithm)
    param.L = L;
    param.infer.symbolMethod = 'pgas';
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_newsticks(data,samples,hyper,param);
    
    % Check that the number of current chains does not exceed maxM
    if(size(samples.seq,1)>param.pgas.maxM)
        % If exceeded, increase the size of param.pgas.particles
        param.pgas.maxM = size(samples.seq,1);
        param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');
    end
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq samples.nest out] = sample_post_Z(data,samples,hyper,param);
    
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
    
    
    %% Genie-aided methods
    if(flag_genieFFBS)
        % If computationally feasible, run FFBS
        param.infer.symbolMethod = 'ffbs';
        param.L = Ltrue;
        [samples_genieFFBS.Z samples_genieFFBS.seq samples_genieFFBS.nest out] = sample_post_Z(data,samples_genieFFBS,hyper,param);
    end
    
    % Run PGAS
    param.infer.symbolMethod = 'pgas';
    param.L = Ltrue;
    [samples_geniePGAS.Z samples_geniePGAS.seq samples_geniePGAS.nest out] = sample_post_Z(data,samples_geniePGAS,hyper,param);

    %% Store current sample
    if(it>param.Niter-param.storeIters)
        samplesAll{it-param.Niter+param.storeIters} = samples;
        if(flag_genieFFBS)
            samplesAll_genieFFBS{it-param.Niter+param.storeIters} = samples_genieFFBS;
        end
        samplesAll_geniePGAS{it-param.Niter+param.storeIters} = samples_geniePGAS;
    end
    
    %% Trace of performance (M_EST and LLH)
    
    % (A) BNP with PGAS method
    M_EST(it) = sum(sum(samples.seq~=0,2)>0);  % Estimated number of transmitters
    param.L = L;
    LLH(it) = compute_llh(data,samples,hyper,param);  % Log-likelihood
    
    % (B) Genie-aided methods
    if(flag_genieFFBS)
        param.L = Ltrue;
        LLH_genieFFBS(it) = compute_llh(data,samples_genieFFBS,hyper,param);  % Log-likelihood
    end
    param.L = Ltrue;
    LLH_geniePGAS(it) = compute_llh(data,samples_geniePGAS,hyper,param);      % Log-likelihood
    
    %% Save temporary result file
    if(mod(it,param.saveCycle)==0)
        save([saveFile '/it' num2str(it) '.mat'],'data','init','initGenie','samples','ADER','SER_ALL','SER_ACT','MMSE','LLH','M_EST','samplesAll',...
                                                 'samples_genieFFBS','ADER_genieFFBS','SER_ALL_genieFFBS','SER_ACT_genieFFBS','LLH_genieFFBS','samplesAll_genieFFBS',...
                                                 'samples_geniePGAS','ADER_geniePGAS','SER_ALL_geniePGAS','SER_ACT_geniePGAS','LLH_geniePGAS','samplesAll_geniePGAS');
        % If successfully saved, delete previous temporary file
        if(exist([saveFile '/it' num2str(it-param.saveCycle) '.mat'],'file'))
            delete([saveFile '/it' num2str(it-param.saveCycle) '.mat']);
        end
    end
end

%% Run BCJR
if(flag_genieBCJR)
    param.L = Ltrue;
    auxConstellation = [0 param.constellation];
    samples_genieBCJR = initGenie;
    [samples_genieBCJR.Z qt_red Simb_red] = bcjr_main(data,samples_genieBCJR,hyper,param);
    samples_genieBCJR.seq = zeros(size(samples_genieBCJR.Z));
    for q=1:length(auxConstellation)
        idx = (abs(samples_genieBCJR.Z-auxConstellation(q))<min(abs(param.constellation))/10);
        samples_genieBCJR.seq(idx) = q-1;
    end
end


%% Final evaluation of performance

% (A) BNP with PGAS method
param.L = L;
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

[ADER(param.Niter+1) SER_ALL(param.Niter+1) SER_ACT(param.Niter+1) MMSE(param.Niter+1) vec_ord rot] = compute_error_rates(data,auxSample,hyper,param);
LLH(param.Niter+1) = compute_llh(data,auxSample,hyper,param);
M_EST(param.Niter+1) = sum(sum(auxSample.seq~=0,2)>0);

% (B) Genie-aided FFBS method
param.L = Ltrue;
if(flag_genieFFBS)
    Zaux = zeros(size(samples_genieFFBS.Z,1),param.T,1+length(param.constellation));
    auxSample = samples_genieFFBS;
    for it=1:param.storeIters
        for t=1:param.T
            for m=1:size(samplesAll_genieFFBS{it}.seq,1)
                Zaux(m,t,samplesAll_genieFFBS{it}.seq(m,t)+1) = 1+Zaux(m,t,samplesAll_genieFFBS{it}.seq(m,t)+1);
            end
        end
    end
    [valnul auxIdx] = max(Zaux,[],3);
    auxConstellation = [0 param.constellation];
    auxSample.seq = auxIdx-1;
    auxSample.Z = auxConstellation(auxIdx);

    [ADER_genieFFBS(param.Niter+1) SER_ALL_genieFFBS(param.Niter+1) SER_ACT_genieFFBS(param.Niter+1) valnul vec_ord rot] = compute_error_rates(data,auxSample,hyper,param,0,0);
    LLH_genieFFBS(param.Niter+1) = compute_llh(data,auxSample,hyper,param);
end

% (C) Genie-aided PGAS method
param.L = Ltrue;
Zaux = zeros(size(samples_geniePGAS.Z,1),param.T,1+length(param.constellation));
auxSample = samples_geniePGAS;
for it=1:param.storeIters
    for t=1:param.T
        for m=1:size(samplesAll_geniePGAS{it}.seq,1)
            Zaux(m,t,samplesAll_geniePGAS{it}.seq(m,t)+1) = 1+Zaux(m,t,samplesAll_geniePGAS{it}.seq(m,t)+1);
        end
    end
end
[valnul auxIdx] = max(Zaux,[],3);
auxConstellation = [0 param.constellation];
auxSample.seq = auxIdx-1;
auxSample.Z = auxConstellation(auxIdx);

[ADER_geniePGAS(param.Niter+1) SER_ALL_geniePGAS(param.Niter+1) SER_ACT_geniePGAS(param.Niter+1) valnul vec_ord rot] = compute_error_rates(data,auxSample,hyper,param,0,0);
LLH_geniePGAS(param.Niter+1) = compute_llh(data,auxSample,hyper,param);

% (D) Genie-aided BCJR method
param.L = Ltrue;
[ADER_genieBCJR(param.Niter+1) SER_ALL_genieBCJR(param.Niter+1) SER_ACT_genieBCJR(param.Niter+1) valnul vec_ord rot] = compute_error_rates(data,samples_genieBCJR,hyper,param,0,0);
LLH_genieBCJR(param.Niter+1) = compute_llh(data,samples_genieBCJR,hyper,param);

%% Save result
save([saveFile '.mat'],'data','init','samples','ADER','SER_ALL','SER_ACT','MMSE','LLH','M_EST','samplesAll',...
                       'samples_genieFFBS','ADER_genieFFBS','SER_ALL_genieFFBS','SER_ACT_genieFFBS','LLH_genieFFBS','samplesAll_genieFFBS',...
                       'samples_geniePGAS','ADER_geniePGAS','SER_ALL_geniePGAS','SER_ACT_geniePGAS','LLH_geniePGAS','samplesAll_geniePGAS',...
                       'ADER_genieBCJR','SER_ALL_genieBCJR','SER_ACT_genieBCJR','LLH_genieBCJR');

%% If successfully saved, detele previous temporary file
if(exist([saveFile '/it' num2str(param.Niter) '.mat'],'file'))
    delete([saveFile '/it' num2str(param.Niter) '.mat']);
end

