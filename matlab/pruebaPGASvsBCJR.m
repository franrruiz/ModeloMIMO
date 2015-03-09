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
param.gen.L_true = [2 2];  % Channel memory for each transmitter (including tap 0)
param.gen.symbol0 = 0;   % Symbols transmitted before (and after) the transmission
param.gen.s2n = 4;                    % Variance of the noise (to generate data)
param.gen.varH = 1;                   % Variance of the channel (to generate data)
param.gen.burstLength = 30*ones(1,param.gen.Nt);  % Mean length of bursts (to generate data)
param.gen.burstLengthStdFactor = 10;
param.gen.sparsityH = 0;

param.Nr = 3;                        % Number of receivers
param.T  = 50;                        % Length of the sequence
M = 4;                                % M for the M-QAM constellation
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = 2;        % Channel memory to be considered during inference

%% Configuration parameters for BCJR and PGAS
param.bcjr.p1 = 0.9;
param.bcjr.p2 = 0.1;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 2000;
param.pgas.returnNsamples = round(0.5*param.pgas.Niter);

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
[Sest qt_red Simb_red] = bcjr_main(data,samples,[],param);
X_PGAS = pgas_main_C(data,samples,[],param);

%% Convert BCJR results into marginals
P_BCJR = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
for t=1:param.T
    for n=1:param.gen.Nt
        a = 0;
        if(param.flag0)
            idx = find(Simb_red(:,n)==0);
            P_BCJR(n,t,1) = sum(qt_red(idx,t));
            a = 1;
        end
        for q=1:length(param.constellation)
            idx = find(Simb_red(:,n)==param.constellation(q));
            P_BCJR(n,t,q+a) = sum(qt_red(idx,t));
        end
    end
end

%% Convert PGAS results into probabilities
M = size(X_PGAS,3);
P_PGAS = zeros(param.gen.Nt,param.T,length(param.constellation)+param.flag0);
for t=1:param.T
    for n=1:param.gen.Nt
        a = 0;
        if(param.flag0)
            P_PGAS(n,t,1) = sum(X_PGAS(n,t,:)==0)/M;
            a = 1;
        end
        for q=1:length(param.constellation)
            P_PGAS(n,t,q+a) = sum(X_PGAS(n,t,:)==param.constellation(q))/M;
        end
    end
end

%% Compute discrepancy:
disp(['Mean deviation: ' num2str(sum(sum(sum(abs(P_BCJR-P_PGAS))))/numel(P_BCJR))]);
[valmax idx] = max(abs(P_BCJR(:)-P_PGAS(:)));
disp(['Max absolute deviation: ' num2str(valmax) ' on a probability of ' num2str(P_BCJR(idx))]);
[valmax idx] = max(abs(P_BCJR(:)-P_PGAS(:))./P_BCJR(:));
disp(['Max relative deviation: ' num2str(valmax)]);

