clear all;
close all;

hohChar = 'C';
L_vec = [1:2 5];
T = 2000;
Nt = 6;
Npart = 3000;
M = 2;
if(hohChar=='B')
    Nr = 11;
    simId = 2;
elseif(hohChar=='C')
    Nr = 12;
    simId = 3;
end

%% Build param
param.constellation = qammod(0:2^M-1,2^M,[],'gray');
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.Nr = Nr;
param.T = T;
param.storeIters = 2000;
param.infer.sampleChannel = 1;

%% Build hyper
hyper.s2h = 0.01;      % E[s2H(r)]=s2h*exp(-lambda*(r-1))
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

%% Storing variables
ALL_MEST = zeros(2,length(L_vec));
ALL_LLH = zeros(2,length(L_vec));
ALL_MMSE = zeros(2,length(L_vec));
ALL_MMSE1 = zeros(2,length(L_vec));

%% Load result files
c = 0;
for L=L_vec
    c = c+1;
    param.L = L;
    
    
    %% Load PGAS
    disp(['Reading PGAS for L=' num2str(L)]);
    load(['/export/clusterdata/franrruiz87/ModeloMIMO/results/wise/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_L' num2str(L) '_lHead0_onOff0_Npart' num2str(Npart) '/itCluster1.mat'],'samples','data','samplesAll');
    
    % Final evaluation of performance
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

    [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samples,hyper,param);
    MMSE1_indiv = zeros(Nt,1);
    for nt=1:Nt
        ll = desp(nt);
        if(ll==0)
            MMSE1_indiv(nt) = sum(abs(data.channel(:,nt,1)-auxSample.H(:,vec_ord(nt),1)/rot(nt)).^2)/Nr;
        elseif(ll>0)
            MMSE1_indiv(nt) = sum(abs(data.channel(:,nt,1)-auxSample.H(:,vec_ord(nt),1+ll)/rot(nt)).^2)/Nr;
        elseif(ll<0)
            error('Why??');
        end
    end
    LLH = compute_llh(data,auxSample,hyper,param);
    M_EST = sum(sum(auxSample.seq~=0,2)>0);
    
    ALL_MEST(1,c) = M_EST;
    ALL_LLH(1,c) = LLH;
    ALL_MMSE(1,c) = mean(MMSE_indiv(1:Nt));
    ALL_MMSE1(1,c) = mean(MMSE1_indiv);
    
    
    %% Load FFBS
    disp(['Reading FFBS for L=' num2str(L)]);
    load(['/export/clusterdata/franrruiz87/ModeloMIMO/results/wise/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_L' num2str(L) '_lHead0_onOff0_FFBS/itCluster1.mat'],'samples','data','samplesAll');
    [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samples,hyper,param);
    
    % Final evaluation of performance
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

    [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samples,hyper,param);
    MMSE1_indiv = zeros(Nt,1);
    for nt=1:Nt
        ll = desp(nt);
        if(ll==0)
            MMSE1_indiv(nt) = sum(abs(data.channel(:,nt,1)-auxSample.H(:,vec_ord(nt),1)/rot(nt)).^2)/Nr;
        elseif(ll>0)
            MMSE1_indiv(nt) = sum(abs(data.channel(:,nt,1)-auxSample.H(:,vec_ord(nt),1+ll)/rot(nt)).^2)/Nr;
        elseif(ll<0)
            error('Why??');
        end
    end
    LLH = compute_llh(data,auxSample,hyper,param);
    M_EST = sum(sum(auxSample.seq~=0,2)>0);
    
    ALL_MEST(2,c) = M_EST;
    ALL_LLH(2,c) = LLH;
    ALL_MMSE(2,c) = mean(MMSE_indiv(1:Nt));
    ALL_MMSE1(2,c) = mean(MMSE1_indiv);    
    
end
