function test_PGAS_C(itCluster)

addpath(genpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab'));
addpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab/auxFunc/mtimesx_20110223');
addpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab/auxFunc/mtimesx_20110223/compiled');

%% Initialize the seed
randn('seed',round(sum(1e5*clock)+itCluster));
rand('seed',round(sum(1e5*clock)+itCluster));

%% Configuration parameters + Parameters of the ground truth (used to generate data)
param.gen.Nt = 2;   % Number of transmitters (to generate data)
param.gen.L_true = [5 5];  % Channel memory for each transmitter (including tap 0)
param.gen.symbol0 = 0;   % Symbols transmitted before (and after) the transmission
param.gen.s2n = 4;                    % Variance of the noise (to generate data)
param.gen.varH = 1;                   % Variance of the channel (to generate data)
param.gen.burstLength = 25*ones(1,param.gen.Nt);  % Mean length of bursts (to generate data)
param.gen.burstLengthStdFactor = 5;
param.gen.sparsityH = 0;

param.Nr = 4;                        % Number of receivers
param.T  = 30;                        % Length of the sequence
M = 4;                                % M for the M-QAM constellation
param.constellation = qammod(0:M-1,M);
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.flag0 = 1;    % Consider symbol 0 as part of the constellation (if false, transmitters are always active)
param.L = 5;        % Channel memory to be considered during inference

%% Configuration parameters for BCJR and PGAS
param.bcjr.p1 = 0.9;
param.bcjr.p2 = 0.1;
param.pgas.N_PF = 3000;
param.pgas.N_PG = 3000;
param.pgas.Niter = 2000;
param.pgas.returnNsamples = round(0.5*param.pgas.Niter);

%% Load data (variable 'data')
load('/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/pruebas/data.mat');

%% Load initialization point (variable 'samples')
load(['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/pruebas/init' num2str(itCluster) '.mat']);

%% Inference
X_PGAS1 = pgas_main(data,samples,[],param);    % This runs matlab PGAS
X_PGAS2 = pgas_main_C(data,samples,[],param);  % This runs C PGAS

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

%% Save results
save(['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/pruebas/results' num2str(itCluster) '.mat'],'X_PGAS1','X_PGAS2','P_PGAS1','P_PGAS2');

