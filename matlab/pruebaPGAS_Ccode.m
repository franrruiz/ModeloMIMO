clear all;
close all;
clc;

addpath bcjr
addpath pgas
addpath pgas/compiled
addpath pgas/mex
addpath sampleFunc
addpath mimoFunc
addpath auxFunc

%% Configuration parameters + Parameters of the ground truth (used to generate data)
param.gen.Nt = 2;   % Number of transmitters (to generate data)
param.gen.L_true = [4 4 2 3 5];  % Channel memory for each transmitter (including tap 0)
param.gen.symbol0 = 0;   % Symbols transmitted before (and after) the transmission
param.gen.s2n = 5;                    % Variance of the noise (to generate data)
param.gen.varH = 1;                   % Variance of the channel (to generate data)
param.gen.burstLength = 10*ones(1,param.gen.Nt);  % Mean length of bursts (to generate data)
param.gen.burstLengthStdFactor = 5;
param.gen.sparsityH = 0;

param.L = 5;        % Channel memory to be considered during inference
param.Nr = 4;                        % Number of receivers
param.T  = 30;                        % Length of the sequence
M = 4;                                % M for the M-QAM constellation
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.header = [];
param.onOffModel = 0;

%% Configuration parameters for BCJR and PGAS
param.bcjr.p1 = 0.9;
param.bcjr.p2 = 0.1;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 2000;
param.pgas.returnNsamples = round(0.5*param.pgas.Niter);
param.pgas.maxM = param.gen.Nt;
param.pgas.particles = zeros(param.pgas.maxM,max(param.pgas.N_PF,param.pgas.N_PG),param.T,'int16');

%% Generate data
data = generate_data_bursts(param);

%% Create struct samples
samples.H = data.channel;
samples.s2y = param.gen.s2n;
samples.am = param.bcjr.p1*ones(param.gen.Nt,1);
samples.bm = param.bcjr.p2*ones(param.gen.Nt,1);
if(param.flag0)
    samples.seq = randint(param.gen.Nt,param.T,[0 M]);
    auxConstellation = [0 param.constellation];
    samples.Z = auxConstellation(samples.seq+1);
end

%% Initialization and inference
tic
X_PGAS1 = pgas_main_matlab(data,samples,[],param);
toc
tic
X_PGAS2 = pgas_main(data,samples,[],param);
toc

%% Convert PGAS results into probabilities
M = size(X_PGAS1,3);
P_PGAS1 = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
P_PGAS2 = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
for t=1:param.T
    for n=1:param.gen.Nt
        a = 0;
        if(param.flag0)
            P_PGAS1(n,t,1) = sum(X_PGAS1(n,t,:)==0)/M;
            P_PGAS2(n,t,1) = sum(X_PGAS2(n,t,:)==0)/M;
            a = 1;
        end
        for q=1:length(param.constellation)
            P_PGAS1(n,t,q+a) = sum(X_PGAS1(n,t,:)==param.constellation(q))/M;
            P_PGAS2(n,t,q+a) = sum(X_PGAS2(n,t,:)==param.constellation(q))/M;
        end
    end
end

%% Compute discrepancy:
disp(['Mean deviation: ' num2str(sum(sum(sum(abs(P_PGAS1-P_PGAS2))))/numel(P_PGAS1))]);
[valmax idx] = max(abs(P_PGAS1(:)-P_PGAS2(:)));
disp(['Max absolute deviation: ' num2str(valmax) ' on a probability of ' num2str(P_PGAS1(idx))]);
[valmax idx] = max(abs(P_PGAS1(:)-P_PGAS2(:))./P_PGAS1(:));
disp(['Max relative deviation: ' num2str(valmax)]);

