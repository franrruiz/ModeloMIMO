clear all;
close all;

plotToFile = 0;

hohChar = 'C';
L_vec = 1:5;
T = 2000;
Nt = 6;
Npart = 3000;
M = 2;
Niter = 30000;
if(hohChar=='B')
    Nr = 11;
    simId = 2;
elseif(hohChar=='C')
    Nr = 12;
    simId = 3;
end

addpath(genpath('/export/clusterdata/franrruiz87/ModeloMIMO/matlab'));

%% Build param
param.constellation = qammod(0:2^M-1,2^M,[],'gray');
param.constellation = param.constellation/sqrt(mean(abs(param.constellation.^2)));
param.Nr = Nr;
param.T = T;
param.storeIters = 2000;
param.infer.sampleChannel = 1;
param.Niter = Niter;

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
BOX_MEST = zeros(2,length(L_vec),param.storeIters);

%% Load result files
c = 0;
for L=L_vec
    c = c+1;
    param.L = L;
    
    
    %% Load PGAS
    if(exist(['/export/clusterdata/franrruiz87/ModeloMIMO/results/wise/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_L' num2str(L) '_lHead0_onOff0_Npart' num2str(Npart) '/itCluster1.mat'],'file'))
        disp(['Reading PGAS for L=' num2str(L)]);
        load(['/export/clusterdata/franrruiz87/ModeloMIMO/results/wise/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_L' num2str(L) '_lHead0_onOff0_Npart' num2str(Npart) '/itCluster1.mat'],'M_EST','samples','data','samplesAll');
    
        % Final evaluation of performance
        LLH = 0;
        MMSE1_indiv = zeros(Nt,1);
        MMSE_tot = 0;
        for it=1:param.storeIters
            [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samplesAll{it},hyper,param);
            for nt=1:Nt
                ll = desp(nt);
                if(ll==0)
                    MMSE1_indiv(nt) = MMSE1_indiv(nt)+sum(abs(data.channel(:,nt,1)-samplesAll{it}.H(:,vec_ord(nt),1)/rot(nt)).^2)/Nr/param.storeIters;
                elseif(ll>0)
                    MMSE1_indiv(nt) = MMSE1_indiv(nt)+sum(abs(data.channel(:,nt,1)-samplesAll{it}.H(:,vec_ord(nt),1+ll)/rot(nt)).^2)/Nr/param.storeIters;
                elseif(ll<0)
                    error('Why??');
                end
            end
            LLH = LLH+compute_llh(data,samplesAll{it},hyper,param)/param.storeIters;
            MMSE_tot = MMSE_tot+mean(MMSE_indiv(1:Nt))/param.storeIters;
        end   
        ALL_MEST(1,c) = mean(M_EST(param.Niter-param.storeIters+1:param.Niter));
        ALL_LLH(1,c) = LLH;
        ALL_MMSE(1,c) = MMSE_tot;
        ALL_MMSE1(1,c) = mean(MMSE1_indiv);
        BOX_MEST(1,c,:) = M_EST(param.Niter-param.storeIters+1:param.Niter);
    end
    
    %% Load FFBS
    if(exist(['/export/clusterdata/franrruiz87/ModeloMIMO/results/wise/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_L' num2str(L) '_lHead0_onOff0_FFBS/itCluster1.mat'],'file'))
        disp(['Reading FFBS for L=' num2str(L)]);
        load(['/export/clusterdata/franrruiz87/ModeloMIMO/results/wise/' num2str(simId) '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_L' num2str(L) '_lHead0_onOff0_FFBS/itCluster1.mat'],'M_EST','samples','data','samplesAll');
        [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samples,hyper,param);

        % Final evaluation of performance
        LLH = 0;
        MMSE1_indiv = zeros(Nt,1);
        MMSE_tot = 0;
        for it=1:param.storeIters
            [ADER SER_ALL SER_ACT MMSE vec_ord rot desp ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samplesAll{it},hyper,param);
            for nt=1:Nt
                ll = desp(nt);
                if(ll==0)
                    MMSE1_indiv(nt) = MMSE1_indiv(nt)+sum(abs(data.channel(:,nt,1)-samplesAll{it}.H(:,vec_ord(nt),1)/rot(nt)).^2)/Nr/param.storeIters;
                elseif(ll>0)
                    MMSE1_indiv(nt) = MMSE1_indiv(nt)+sum(abs(data.channel(:,nt,1)-samplesAll{it}.H(:,vec_ord(nt),1+ll)/rot(nt)).^2)/Nr/param.storeIters;
                elseif(ll<0)
                    error('Why??');
                end
            end
            LLH = LLH+compute_llh(data,samplesAll{it},hyper,param)/param.storeIters;
            MMSE_tot = MMSE_tot+mean(MMSE_indiv(1:Nt))/param.storeIters;
        end   
        ALL_MEST(2,c) = mean(M_EST(param.Niter-param.storeIters+1:param.Niter));
        ALL_LLH(2,c) = LLH;
        ALL_MMSE(2,c) = MMSE_tot;
        ALL_MMSE1(2,c) = mean(MMSE1_indiv);
        BOX_MEST(2,c,:) = M_EST(param.Niter-param.storeIters+1:param.Niter);
    end
end


%% Boxplot
figure;
bar(L_vec,ALL_MEST');
set(gca,'FontSize',14);
legend('iFDM','iFHMM','Location','NorthEast');
xlabel('L');
ylabel('M_+');
grid on;
if(plotToFile)
    figurapdf(4.5,3);  % Before: 3,2
    print('-dpdf',['CommMestNIPS.pdf']);
end




