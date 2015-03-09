clear all;
close all;
clc;

addpath bcjr
addpath pgas


%% Configuration parameters + Parameters of the ground truth (used to generate data)
param.gen.Nt = 3;   % Number of transmitters (to generate data)
param.gen.L_true = 1*ones(1,param.gen.Nt);  % Channel memory for each transmitter (including tap 0)
param.gen.symbol0 = 0;   % Symbols transmitted before (and after) the transmission
param.gen.s2n = 0.5;                    % Variance of the noise (to generate data)
param.gen.varH = 2;                   % Variance of the channel (to generate data)
param.gen.burstLength = 150*ones(1,param.gen.Nt);  % Mean length of bursts (to generate data)

param.Nr = 25;                        % Number of receivers
param.T  = 500;                       % Length of the sequence
M = 4;                                % M for the M-QAM constellation
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = 1;        % Channel memory to be considered during inference

%% Configuration parameters for BCJR and PGAS
param.bcjr.p1 = 0.95;
param.bcjr.p2 = 0.05;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 1;
param.pgas.returnNsamples = 1;

%% Configuration parameters for BNP
param.bnp.betaSlice1 = 0.1;
param.bnp.betaSlice2 = 5;
param.bnp.maxMnew = 15;

%% Generate data
data = generate_data_bursts(param);

%% Hyperparameters
hyper.s2h = 1;      % Variance of H_r ~ CN(0,s2h*exp(-lambda*r)*I)
hyper.lambda = 0.1; % Variance of H_r ~ CN(0,s2h*exp(-lambda*r)*I)
hyper.alpha = 1;    % Concentration parameter for Z ~ IBP(alpha)
hyper.gamma1 = 1;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.gamma2 = 1;   % Parameter for bm ~ Beta(gamma1,gamma2)
hyper.tau = 1;      % Parameter for s2y ~ IG(tau,nu)
hyper.nu = 1;       % Parameter for s2y ~ IG(tau,nu)

%% Initialization
samples.H = zeros(param.Nr,param.gen.Nt,param.L);
for ll=1:param.L
    samples.H(:,:,ll) = sqrt(hyper.s2h*exp(-hyper.lambda*(ll-1)))*(randn(param.Nr,param.gen.Nt,1)+1i*randn(param.Nr,param.gen.Nt,1));
end
samples.s2y = param.gen.s2n;
samples.am = 0.95*ones(param.gen.Nt,1);
samples.bm = 0.05*ones(param.gen.Nt,1);
samples.Z = [];
samples.seq = [];
samples.slice = 0;

%% Inference
for i=1:5000
    i
    %samples.seq = [];
    
    % Step 1)
    % -Sample the slice variable
    samples.slice = sample_post_slice(data,samples,hyper,param);
    % -Sample new sticks (and the corresponding new parameters)
    samples = sample_post_newsticks(data,samples,hyper,param);
    
    % Step 2)
    % -Sample the symbols Z
    [samples.Z samples.seq] = pgas_main(data,samples,hyper,param);
    
    % Step 3)
    
    % Step 4)
    % -Sample the transition probabilities
    
    
    % Step 5)
    % -Sample the channel H
    samples.H = sample_post_H(data,samples,hyper,param);
    % -Sample the noise variance
    %samples.s2y = sample_post_s2y(data,samples,hyper,param);
end




