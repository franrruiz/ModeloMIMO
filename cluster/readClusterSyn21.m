clear all;
close all;
clc
addpath aboxplot;

maxItCluster = 50;
Niter =  20000;
maxNt = 20;

%% Scenario definition


%%%%%%%%%%%%%%%
plotToFile = 1;
flagPlotGenie = 1;
thrSER = 0.1;
%%%%%%%%%%%%%%%

%Base scenario
T = 1000;
Nt = 5;
SNR = -3;
Nr = 20;
L = 1;
Ltrue = 1;
M = 2;
lHead = 0;
onOffModel = 0;
Nparticles = 3000;
MAX_Nt = 10;

% Sweep and hold variables:
  % L = 1:
  %sweep_var = 'SNR'; sweep_vec = -12:3:0; simId = 21; etiquetaX='-10log(\sigma_y^2)'; lugar='NorthEast'; hold_var='Nr'; hold_vec=20; subcarpeta = 'L1/';
  %sweep_var = 'Nt'; sweep_vec = 2:1:6; simId = 21; etiquetaX='#Transmitters'; lugar='NorthWest'; hold_var='SNR'; hold_vec=-3; subcarpeta = 'L1/';
  sweep_var = 'Nr'; sweep_vec = [2:2:10 15:5:30]; simId = 21; etiquetaX='#Receivers'; lugar='NorthEast'; hold_var='SNR'; hold_vec=-3; subcarpeta = 'L1/';
  %sweep_var = 'M'; sweep_vec = 2:1:7; simId = 21; etiquetaX='log_2|A|'; lugar='SouthEast'; hold_var='SNR'; hold_vec=12; Niter = 100000; subcarpeta = 'L1/'; %%% CUIDADO

  % L = 5:
  %sweep_var = 'SNR'; sweep_vec = -15:3:3; simId = 21; etiquetaX='-10log(\sigma_y^2)'; lugar='NorthEast'; hold_var='Ltrue'; hold_vec=5; L=5; subcarpeta = 'L5/';

  % Sweep Ltrue:
  %sweep_var = 'Ltrue'; sweep_vec = 1:2:7; simId = 21; etiquetaX='L'; lugar='NorthEast'; hold_var='SNR'; hold_vec=-9; subcarpeta = 'sweep_Ltrue/';
  
  % Sweep L:
  %sweep_var = 'L'; sweep_vec = 1:1:3; simId = 21; etiquetaX='L'; lugar='NorthEast'; hold_var='Ltrue'; hold_vec=1; simId = 23; subcarpeta = 'sweep_L/';

marcadores = {'+','^','o','x'};
estilos = {'-','--',':','-.'};
colores = {'m','b','k','r'};
leyenda = {'IFFSM','PGAS','BCJR','FFBS'};
anchoSm = 4; altoSm = 3;  % Before: 3,2
anchoLar = 5; altoLar = 3.5;

%% Initialize metrics of interest
ALL_ADER = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MMSE = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_LLH = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MEST = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_RECOVERED = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));

ALL_ADER_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_MMSE_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));

ALL_ADER_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_MMSE_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_ADER_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_MMSE_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_ADER_FFBS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_FFBS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_FFBS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_MMSE_FFBS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));

ALL_ADER_PGAS_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_PGAS_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_PGAS_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_ADER_BCJR_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_BCJR_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_BCJR_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_ADER_FFBS_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_FFBS_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_FFBS_indiv = -ones(maxItCluster,maxNt,length(sweep_vec),length(hold_vec));

%% Load results
idxItClusterNotFound = [];
for sweepV=sweep_vec
    idxSweep = find(sweep_vec==sweepV);
    
    % Go back to base configuration
    Taux = T;
    Ntaux = Nt;
    SNRaux = SNR;
    Nraux = Nr;
    Laux = L;
    Ltrueaux = Ltrue;
    Maux = M;
    lHeadaux = lHead;
    onOffModelaux = onOffModel;
    
    % Modify the sweep variable
    if(strcmp(sweep_var,'T'))
        Taux = sweepV;
    elseif(strcmp(sweep_var,'Nt'))
        Ntaux = sweepV;
    elseif(strcmp(sweep_var,'SNR'))
        SNRaux = sweepV;
    elseif(strcmp(sweep_var,'Nr'))
        Nraux = sweepV;
    elseif(strcmp(sweep_var,'L'))
        Laux = sweepV;
    elseif(strcmp(sweep_var,'Ltrue'))
        Ltrueaux = sweepV;
        Laux = sweepV;
    elseif(strcmp(sweep_var,'M'))
        Maux = sweepV;
    elseif(strcmp(sweep_var,'lHead'))
        lHeadaux = sweepV;
    elseif(strcmp(sweep_var,'onOffModel'))
        onOffModelaux = sweepV;
    else
        error('Not valid value for sweep_var')
    end
    
    for holdV=hold_vec
        idxHold = find(hold_vec==holdV);
        
        % Modify the hold variable
        if(strcmp(hold_var,'T'))
            Taux = holdV;
        elseif(strcmp(hold_var,'Nt'))
            Ntaux = holdV;
        elseif(strcmp(hold_var,'SNR'))
            SNRaux = holdV;
        elseif(strcmp(hold_var,'Nr'))
            Nraux = holdV;
        elseif(strcmp(hold_var,'L'))
            Laux = holdV;
        elseif(strcmp(hold_var,'Ltrue'))
            Ltrueaux = holdV;
        elseif(strcmp(hold_var,'M'))
            Maux = holdV;
        elseif(strcmp(hold_var,'lHead'))
            lHeadaux = holdV;
        elseif(strcmp(hold_var,'onOffModel'))
            onOffModelaux = holdV;
        else
            error('Not valid value for hold_var')
        end
        
        % Load all itCluster's for the given configuration
        for itCluster=1:maxItCluster
            
            
            saveFile = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId)...
                        '/T' num2str(Taux) '_Nt' num2str(Ntaux) '_Nr' num2str(Nraux) '_M' num2str(Maux) '_Ltrue' num2str(Ltrueaux) '_L' num2str(Laux) '_SNR' num2str(SNRaux) '_lHead' num2str(lHeadaux), '_onOff' num2str(onOffModelaux) '_Npart' num2str(Nparticles) ...
                        '/itCluster' num2str(itCluster) '.mat'];
                    
            if(exist(saveFile,'file'))
                clear ADER_PGAS4 ADER_BCJR4 ADER_FFBS4
                load(saveFile,'ADER','MMSE','SER_ACT','SER_ALL','LLH','M_EST','*_indiv','*_PGAS4','*_BCJR4','*_FFBS4');
                
                ALL_ADER(itCluster,:,idxSweep,idxHold) = ADER;
                ALL_MMSE(itCluster,:,idxSweep,idxHold) = MMSE;
                ALL_LLH(itCluster,:,idxSweep,idxHold) = LLH;
                ALL_SER_ALL(itCluster,:,idxSweep,idxHold) = SER_ALL;
                ALL_SER_ACT(itCluster,:,idxSweep,idxHold) = SER_ACT;
                ALL_MEST(itCluster,:,idxSweep,idxHold) = M_EST;
                ALL_RECOVERED(itCluster,1,idxSweep,idxHold) = sum(SER_ALL_indiv<thrSER);
                
                ALL_ADER_indiv(itCluster,1:length(ADER_indiv),idxSweep,idxHold) = ADER_indiv;
                ALL_MMSE_indiv(itCluster,1:length(MMSE_indiv),idxSweep,idxHold) = MMSE_indiv;
                ALL_SER_ALL_indiv(itCluster,1:length(SER_ALL_indiv),idxSweep,idxHold) = SER_ALL_indiv;
                ALL_SER_ACT_indiv(itCluster,1:length(SER_ACT_indiv),idxSweep,idxHold) = SER_ACT_indiv;
                
                if(exist('ADER_PGAS4','var'))
                    ALL_ADER_PGAS(itCluster,1,idxSweep,idxHold) = ADER_PGAS4;
                    ALL_SER_ALL_PGAS(itCluster,1,idxSweep,idxHold) = SER_ALL_PGAS4;
                    ALL_SER_ACT_PGAS(itCluster,1,idxSweep,idxHold) = SER_ACT_PGAS4;
                    ALL_MMSE_PGAS(itCluster,1,idxSweep,idxHold) = MMSE_PGAS4;
                    
                    ALL_ADER_PGAS_indiv(itCluster,1:length(ADER_indiv),idxSweep,idxHold) = ADER_PGAS4_indiv;
                    ALL_SER_ALL_PGAS_indiv(itCluster,1:length(SER_ALL_indiv),idxSweep,idxHold) = SER_ALL_PGAS4_indiv;
                    ALL_SER_ACT_PGAS_indiv(itCluster,1:length(SER_ACT_indiv),idxSweep,idxHold) = SER_ACT_PGAS4_indiv;
                end
                if(exist('ADER_BCJR4','var'))
                    ALL_ADER_BCJR(itCluster,1,idxSweep,idxHold) = ADER_BCJR4;
                    ALL_SER_ALL_BCJR(itCluster,1,idxSweep,idxHold) = SER_ALL_BCJR4;
                    ALL_SER_ACT_BCJR(itCluster,1,idxSweep,idxHold) = SER_ACT_BCJR4;
                    ALL_MMSE_BCJR(itCluster,1,idxSweep,idxHold) = MMSE_BCJR4;
                    
                    ALL_ADER_BCJR_indiv(itCluster,1:length(ADER_indiv),idxSweep,idxHold) = ADER_BCJR4_indiv;
                    ALL_SER_ALL_BCJR_indiv(itCluster,1:length(SER_ALL_indiv),idxSweep,idxHold) = SER_ALL_BCJR4_indiv;
                    ALL_SER_ACT_BCJR_indiv(itCluster,1:length(SER_ACT_indiv),idxSweep,idxHold) = SER_ACT_BCJR4_indiv;
                end
                if(exist('ADER_FFBS4','var'))
                    ALL_ADER_FFBS(itCluster,1,idxSweep,idxHold) = ADER_FFBS4;
                    ALL_SER_ALL_FFBS(itCluster,1,idxSweep,idxHold) = SER_ALL_FFBS4;
                    ALL_SER_ACT_FFBS(itCluster,1,idxSweep,idxHold) = SER_ACT_FFBS4;
                    ALL_MMSE_FFBS(itCluster,1,idxSweep,idxHold) = MMSE_FFBS4;
                    
                    ALL_ADER_FFBS_indiv(itCluster,1:length(ADER_indiv),idxSweep,idxHold) = ADER_FFBS4_indiv;
                    ALL_SER_ALL_FFBS_indiv(itCluster,1:length(SER_ALL_indiv),idxSweep,idxHold) = SER_ALL_FFBS4_indiv;
                    ALL_SER_ACT_FFBS_indiv(itCluster,1:length(SER_ACT_indiv),idxSweep,idxHold) = SER_ACT_FFBS4_indiv;
                end
            else
                idxItClusterNotFound = unique([idxItClusterNotFound itCluster]);
                disp(['File ' saveFile ' not found']);
            end
        end
        
    end
end

%% Remove not found itCluster's
%idxItClusterNotFound = unique([idxItClusterNotFound 19]);
if(~isempty(idxItClusterNotFound))
    disp(['WARNING: Removing ' num2str(length(idxItClusterNotFound)) ' simulations. Total is now: ' num2str(maxItCluster-length(idxItClusterNotFound))]);
end
ALL_ADER(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE(idxItClusterNotFound,:,:,:) = [];
ALL_LLH(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT(idxItClusterNotFound,:,:,:) = [];
ALL_MEST(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_PGAS(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_PGAS(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_PGAS(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE_PGAS(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_BCJR(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_BCJR(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_BCJR(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE_BCJR(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE_FFBS(idxItClusterNotFound,:,:,:) = [];

ALL_RECOVERED(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_indiv(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_PGAS_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_PGAS_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_PGAS_indiv(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_BCJR_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_BCJR_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_BCJR_indiv(idxItClusterNotFound,:,:,:) = [];

ALL_ADER_FFBS_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_FFBS_indiv(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_FFBS_indiv(idxItClusterNotFound,:,:,:) = [];

%% Build idxGood_M, idxBad_M
Ntaux = Nt;
if(strcmp(sweep_var,'Nt'))
    vecAux = zeros(1,1,length(sweep_vec));
    vecAux(:) = sweep_vec;
    idxGood_M = (ALL_MEST(:,end,:,:)==repmat(vecAux,[maxItCluster-length(idxItClusterNotFound),1,1,length(hold_vec)]));
    idxBad_M = (ALL_MEST(:,end,:,:)==repmat(vecAux,[maxItCluster-length(idxItClusterNotFound),1,1,length(hold_vec)]));
elseif(strcmp(hold_var,'Nt'))
    vecAux = zeros(1,1,1,length(hold_vec));
    vecAux(:) = hold_vec;
    idxGood_M = (ALL_MEST(:,end,:,:)==repmat(vecAux,[maxItCluster-length(idxItClusterNotFound),1,length(sweep_vec),1]));
    idxBad_M = (ALL_MEST(:,end,:,:)==repmat(vecAux,[maxItCluster-length(idxItClusterNotFound),1,length(sweep_vec),1]));
else
    idxGood_M = ((ALL_MEST(:,end,:,:)==repmat(Ntaux,[maxItCluster-length(idxItClusterNotFound),1,length(sweep_vec),length(hold_vec)])));
    idxBad_M = ((ALL_MEST(:,end,:,:)==repmat(Ntaux,[maxItCluster-length(idxItClusterNotFound),1,length(sweep_vec),length(hold_vec)])));
end

%% For fixed values of the hold variable, plot the metrics vs the sweep variable
for holdV=hold_vec
    idxHold = find(hold_vec==holdV);
    
    % Go back to base configuration
    Taux = T;
    Ntaux = Nt;
    SNRaux = SNR;
    Nraux = Nr;
    Laux = L;
    Ltrueaux = Ltrue;
    Maux = M;
    lHeadaux = lHead;
    onOffModelaux = onOffModel;
    
    if(strcmp(hold_var,'Nt'))
        Ntaux = holdV;
    elseif(strcmp(hold_var,'Ltrue'))
        Ltrueaux = holdV;
    elseif(strcmp(hold_var,'L'))
        Laux = holdV;
    end
    
    % Load the results of interest
    avgADER = zeros(1,length(sweep_vec));
    avgSER_ACT = zeros(1,length(sweep_vec));
    avgSER_ALL = zeros(1,length(sweep_vec));
    avgMSE = zeros(1,length(sweep_vec));
    avgMEST = zeros(1,length(sweep_vec));
    boxMEST = zeros(maxItCluster-length(idxItClusterNotFound),length(sweep_vec));
    stdADER = zeros(1,length(sweep_vec));
    stdSER_ACT = zeros(1,length(sweep_vec));
    stdSER_ALL = zeros(1,length(sweep_vec));
    stdMSE = zeros(1,length(sweep_vec));
    stdMEST = zeros(1,length(sweep_vec));
    
    avgADER_PGAS = zeros(1,length(sweep_vec));
    avgSER_ACT_PGAS = zeros(1,length(sweep_vec));
    avgSER_ALL_PGAS = zeros(1,length(sweep_vec));
    stdADER_PGAS = zeros(1,length(sweep_vec));
    stdSER_ACT_PGAS = zeros(1,length(sweep_vec));
    stdSER_ALL_PGAS = zeros(1,length(sweep_vec));

    avgADER_BCJR = zeros(1,length(sweep_vec));
    avgSER_ACT_BCJR = zeros(1,length(sweep_vec));
    avgSER_ALL_BCJR = zeros(1,length(sweep_vec));
    stdADER_BCJR = zeros(1,length(sweep_vec));
    stdSER_ACT_BCJR = zeros(1,length(sweep_vec));
    stdSER_ALL_BCJR = zeros(1,length(sweep_vec));
    
    avgADER_FFBS = zeros(1,length(sweep_vec));
    avgSER_ACT_FFBS = zeros(1,length(sweep_vec));
    avgSER_ALL_FFBS = zeros(1,length(sweep_vec));
    stdADER_FFBS = zeros(1,length(sweep_vec));
    stdSER_ACT_FFBS = zeros(1,length(sweep_vec));
    stdSER_ALL_FFBS = zeros(1,length(sweep_vec));
    
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Ntaux = sweepV;
        elseif(strcmp(sweep_var,'Ltrue'))
            Ltrueaux = sweepV;
            Laux = sweepV;
        end
        
        auxADER = [];
        auxSER_ALL = [];
        auxSER_ACT = [];
        auxMMSE = [];
        auxADER_PGAS = [];
        auxSER_ALL_PGAS = [];
        auxSER_ACT_PGAS = [];
        auxADER_BCJR = [];
        auxSER_ALL_BCJR = [];
        auxSER_ACT_BCJR = [];
        auxADER_FFBS = [];
        auxSER_ALL_FFBS = [];
        auxSER_ACT_FFBS = [];
        for ii=1:maxItCluster-length(idxItClusterNotFound)
            idxGood_indiv = (ALL_SER_ALL_indiv(ii,:,idxSweep,idxHold)<thrSER & ALL_SER_ALL_indiv(ii,:,idxSweep,idxHold)>=0);
            auxADER = [auxADER ALL_ADER_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ALL = [auxSER_ALL ALL_SER_ALL_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ACT = [auxSER_ACT ALL_SER_ACT_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxMMSE = [auxMMSE ALL_MMSE_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            
            idxGood_indiv = (ALL_SER_ALL_PGAS_indiv(ii,:,idxSweep,idxHold)<thrSER & ALL_SER_ALL_PGAS_indiv(ii,:,idxSweep,idxHold)>=0);
            auxADER_PGAS = [auxADER_PGAS ALL_ADER_PGAS_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ALL_PGAS = [auxSER_ALL_PGAS ALL_SER_ALL_PGAS_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ACT_PGAS = [auxSER_ACT_PGAS ALL_SER_ACT_PGAS_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            
            idxGood_indiv = (ALL_SER_ALL_BCJR_indiv(ii,:,idxSweep,idxHold)<thrSER & ALL_SER_ALL_BCJR_indiv(ii,:,idxSweep,idxHold)>=0);
            auxADER_BCJR = [auxADER_BCJR ALL_ADER_BCJR_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ALL_BCJR = [auxSER_ALL_BCJR ALL_SER_ALL_BCJR_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ACT_BCJR = [auxSER_ACT_BCJR ALL_SER_ACT_BCJR_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            
            idxGood_indiv = (ALL_SER_ALL_FFBS_indiv(ii,:,idxSweep,idxHold)<thrSER & ALL_SER_ALL_FFBS_indiv(ii,:,idxSweep,idxHold)>=0);
            auxADER_FFBS = [auxADER_FFBS ALL_ADER_FFBS_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ALL_FFBS = [auxSER_ALL_FFBS ALL_SER_ALL_FFBS_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
            auxSER_ACT_FFBS = [auxSER_ACT_FFBS ALL_SER_ACT_FFBS_indiv(ii,idxGood_indiv,idxSweep,idxHold)];
        end
        avgADER(idxSweep) = mean(auxADER);
        stdADER(idxSweep) = std(auxADER);
        avgSER_ACT(idxSweep) = mean(auxSER_ACT);
        stdSER_ACT(idxSweep) = std(auxSER_ACT);
        avgSER_ALL(idxSweep) = mean(auxSER_ALL);
        stdSER_ALL(idxSweep) = std(auxSER_ALL);
        avgMEST(idxSweep) = mean(ALL_MEST(:,end,idxSweep,idxHold));
        stdMEST(idxSweep) = std(ALL_MEST(:,end,idxSweep,idxHold));
        avgMSE(idxSweep) = mean(auxMMSE);
        stdMSE(idxSweep) = std(auxMMSE);
        boxMEST(:,idxSweep) = abs(ALL_MEST(:,end,idxSweep,idxHold));
        
        avgADER_PGAS(idxSweep) = mean(auxADER_PGAS);
        stdADER_PGAS(idxSweep) = std(auxADER_PGAS);
        avgSER_ACT_PGAS(idxSweep) = mean(auxSER_ACT_PGAS);
        stdSER_ACT_PGAS(idxSweep) = std(auxSER_ACT_PGAS);
        avgSER_ALL_PGAS(idxSweep) = mean(auxSER_ALL_PGAS);
        stdSER_ALL_PGAS(idxSweep) = std(auxSER_ALL_PGAS);
        
        avgADER_BCJR(idxSweep) = mean(auxADER_BCJR);
        stdADER_BCJR(idxSweep) = std(auxADER_BCJR);
        avgSER_ACT_BCJR(idxSweep) = mean(auxSER_ACT_BCJR);
        stdSER_ACT_BCJR(idxSweep) = std(auxSER_ACT_BCJR);
        avgSER_ALL_BCJR(idxSweep) = mean(auxSER_ALL_BCJR);
        stdSER_ALL_BCJR(idxSweep) = std(auxSER_ALL_BCJR);
        
        avgADER_FFBS(idxSweep) = mean(auxADER_FFBS);
        stdADER_FFBS(idxSweep) = std(auxADER_FFBS);
        avgSER_ACT_FFBS(idxSweep) = mean(auxSER_ACT_FFBS);
        stdSER_ACT_FFBS(idxSweep) = std(auxSER_ACT_FFBS);
        avgSER_ALL_FFBS(idxSweep) = mean(auxSER_ALL_FFBS);
        stdSER_ALL_FFBS(idxSweep) = std(auxSER_ALL_FFBS);
        
%         idxGood = 1:maxItCluster; %find((ALL_MEST(:,end,idxSweep,idxHold)==Ntaux));
%         idxBad = 1:maxItCluster;  %find((ALL_MEST(:,end,idxSweep,idxHold)~=Ntaux));
%         
%         avgADER(idxSweep) = mean(ALL_ADER(idxGood,end,idxSweep,idxHold));
%         stdADER(idxSweep) = std(ALL_ADER(idxGood,end,idxSweep,idxHold));
%         avgSER_ACT(idxSweep) = mean(ALL_SER_ACT(idxGood,end,idxSweep,idxHold));
%         stdSER_ACT(idxSweep) = std(ALL_SER_ACT(idxGood,end,idxSweep,idxHold));
%         avgSER_ALL(idxSweep) = mean(ALL_SER_ALL(idxGood,end,idxSweep,idxHold));
%         stdSER_ALL(idxSweep) = std(ALL_SER_ALL(idxGood,end,idxSweep,idxHold));
%         avgMEST(idxSweep) = mean(ALL_MEST(idxGood,end,idxSweep,idxHold));
%         stdMEST(idxSweep) = std(ALL_MEST(idxGood,end,idxSweep,idxHold));
%         avgMSE(idxSweep) = mean(ALL_MMSE(idxGood,end,idxSweep,idxHold));
%         stdMSE(idxSweep) = std(ALL_MMSE(idxGood,end,idxSweep,idxHold));
%         boxMEST(:,idxSweep) = abs(ALL_MEST(:,end,idxSweep,idxHold));
%         
%         avgADER_PGAS(idxSweep) = mean(ALL_ADER_PGAS(idxGood,1,idxSweep,idxHold));
%         stdADER_PGAS(idxSweep) = std(ALL_ADER_PGAS(idxGood,1,idxSweep,idxHold));
%         avgSER_ACT_PGAS(idxSweep) = mean(ALL_SER_ACT_PGAS(idxGood,1,idxSweep,idxHold));
%         stdSER_ACT_PGAS(idxSweep) = std(ALL_SER_ACT_PGAS(idxGood,1,idxSweep,idxHold));
%         avgSER_ALL_PGAS(idxSweep) = mean(ALL_SER_ALL_PGAS(idxGood,1,idxSweep,idxHold));
%         stdSER_ALL_PGAS(idxSweep) = std(ALL_SER_ALL_PGAS(idxGood,1,idxSweep,idxHold));
%         
%         avgADER_BCJR(idxSweep) = mean(ALL_ADER_BCJR(idxGood,1,idxSweep,idxHold));
%         stdADER_BCJR(idxSweep) = std(ALL_ADER_BCJR(idxGood,1,idxSweep,idxHold));
%         avgSER_ACT_BCJR(idxSweep) = mean(ALL_SER_ACT_BCJR(idxGood,1,idxSweep,idxHold));
%         stdSER_ACT_BCJR(idxSweep) = std(ALL_SER_ACT_BCJR(idxGood,1,idxSweep,idxHold));
%         avgSER_ALL_BCJR(idxSweep) = mean(ALL_SER_ALL_BCJR(idxGood,1,idxSweep,idxHold));
%         stdSER_ALL_BCJR(idxSweep) = std(ALL_SER_ALL_BCJR(idxGood,1,idxSweep,idxHold));
%         
%         avgADER_FFBS(idxSweep) = mean(ALL_ADER_FFBS(idxGood,1,idxSweep,idxHold));
%         stdADER_FFBS(idxSweep) = std(ALL_ADER_FFBS(idxGood,1,idxSweep,idxHold));
%         avgSER_ACT_FFBS(idxSweep) = mean(ALL_SER_ACT_FFBS(idxGood,1,idxSweep,idxHold));
%         stdSER_ACT_FFBS(idxSweep) = std(ALL_SER_ACT_FFBS(idxGood,1,idxSweep,idxHold));
%         avgSER_ALL_FFBS(idxSweep) = mean(ALL_SER_ALL_FFBS(idxGood,1,idxSweep,idxHold));
%         stdSER_ALL_FFBS(idxSweep) = std(ALL_SER_ALL_FFBS(idxGood,1,idxSweep,idxHold));
    end
    
    % Plot ADER
    figure;
    h = semilogy(sweep_vec,avgADER,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('ADER');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    set(gca,'YLim',[10^floor(log10(min(avgADER))) 10^ceil(log10(max(avgADER)))]);
    if(simId==2 && strcmp(sweep_var,'SNR'))
        set(gca,'XTickLabel',sweep_vec-13);
    else
        set(gca,'XTickLabel',sweep_vec);
    end
    grid on;
    hold on;
    h = plot(sweep_vec,avgADER_PGAS,'Marker',marcadores{2},'LineStyle',estilos{2},'Color',colores{2});
    if(~flagPlotGenie)
        set(h,'Visible','off');
    end
    idxBCJR = [];
    idxFFBS = [];
    if(strcmp(sweep_var,'Nt'))
        idxBCJR = find((1+2^Maux).^(2*Laux*sweep_vec)<1e6);
        if((1+2^Maux)^(2*Laux)<1e6)
            idxFFBS = 1:length(sweep_vec);
        end
    elseif(strcmp(sweep_var,'M'))
        idxBCJR = find((1+2.^sweep_vec).^(2*Laux*Ntaux)<1e6);
        idxFFBS = find((1+2.^sweep_vec).^(2*Laux)<1e6);
    elseif(strcmp(sweep_var,'L'))
        idxBCJR = find((1+2^Maux).^(2*sweep_vec*Ntaux)<1e6);
        idxFFBS = find((1+2^Maux).^(2*sweep_vec)<1e6);
    elseif(strcmp(sweep_var,'Ltrue'))
        idxBCJR = find((1+2^Maux).^(2*sweep_vec*Ntaux)<1e6);
        idxFFBS = find((1+2^Maux).^(2*sweep_vec)<1e6);
    else
        if((1+2^Maux)^(2*Laux*Ntaux)<1e6)
            idxBCJR = 1:length(sweep_vec);
        end
        if((1+2^Maux)^(2*Laux)<1e6)
            idxFFBS = 1:length(sweep_vec);
        end
    end
    if(~isempty(idxBCJR))
        hold on;
        h = plot(sweep_vec(idxBCJR),avgADER_BCJR(idxBCJR),'Marker',marcadores{3},'LineStyle',estilos{3},'Color',colores{3});
        if(~flagPlotGenie)
            set(h,'Visible','off');
        end
    end
    if(~isempty(idxFFBS))
        hold on;
        h = plot(sweep_vec(idxFFBS),avgADER_FFBS(idxFFBS),'Marker',marcadores{4},'LineStyle',estilos{4},'Color',colores{4});
        if(~flagPlotGenie)
            set(h,'Visible','off');
        end
    end
    if(flagPlotGenie)
        if(isempty(idxBCJR) && ~isempty(idxFFBS))
            legend(leyenda([1 2 4]),'Location',lugar);
        elseif(~isempty(idxBCJR) && ~isempty(idxFFBS))
            legend(leyenda([1 2 3 4]),'Location',lugar);
        elseif(~isempty(idxBCJR) && isempty(idxFFBS))
            legend(leyenda([1 2 3]),'Location',lugar);
        else
            legend(leyenda([1 2]),'Location',lugar);
        end
    end
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn21/' subcarpeta 'ADER_' sweep_var '_s.pdf']);
    end
    
    % Plot SER_ALL
    figure;
    h = semilogy(sweep_vec,avgSER_ALL,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('SER');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    set(gca,'YLim',[10^floor(log10(min(avgSER_ALL))) 10^ceil(log10(max(avgSER_ALL)))]);
    if(simId==2 && strcmp(sweep_var,'SNR'))
        set(gca,'XTickLabel',sweep_vec-13);
    else
        set(gca,'XTickLabel',sweep_vec);
    end
    hold on;
    grid on;
    h = plot(sweep_vec,avgSER_ALL_PGAS,'Marker',marcadores{2},'LineStyle',estilos{2},'Color',colores{2});
    if(~flagPlotGenie)
        set(h,'Visible','off');
    end
    if(~isempty(idxBCJR))
        hold on;
        h = plot(sweep_vec(idxBCJR),avgSER_ALL_BCJR(idxBCJR),'Marker',marcadores{3},'LineStyle',estilos{3},'Color',colores{3});
        if(~flagPlotGenie)
            set(h,'Visible','off');
        end
    end
    if(~isempty(idxFFBS))
        hold on;
        h = plot(sweep_vec(idxFFBS),avgSER_ALL_FFBS(idxFFBS),'Marker',marcadores{4},'LineStyle',estilos{4},'Color',colores{4});
        if(~flagPlotGenie)
            set(h,'Visible','off');
        end
    end
    if(flagPlotGenie)
        if(isempty(idxBCJR) && ~isempty(idxFFBS))
            legend(leyenda([1 2 4]),'Location',lugar);
        elseif(~isempty(idxBCJR) && ~isempty(idxFFBS))
            legend(leyenda([1 2 3 4]),'Location',lugar);
        elseif(~isempty(idxBCJR) && isempty(idxFFBS))
            legend(leyenda([1 2 3]),'Location',lugar);
        else
            legend(leyenda([1 2]),'Location',lugar);
        end
    end
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn21/' subcarpeta 'SER_' sweep_var '_s.pdf']);
    end
    
    % Plot MSE
    figure;
    h = semilogy(sweep_vec,avgMSE,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('MSE');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    set(gca,'YLim',[10^floor(log10(min(avgMSE))) 10^ceil(log10(max(avgMSE)))]);
    if(simId==2 && strcmp(sweep_var,'SNR'))
        set(gca,'XTickLabel',sweep_vec-13);
    else
        set(gca,'XTickLabel',sweep_vec);
    end
    grid on;
    %legend(leyenda(1));
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn21/' subcarpeta 'MSE_' sweep_var '_s.pdf']);
    end
    
    % Plot DEP
    figure;
    aux_vec = Ntaux*ones(size(sweep_vec));
    if(strcmp(sweep_var,'Nt'))
        aux_vec = sweep_vec;
    end
    h = plot(sweep_vec,1-sum(boxMEST==repmat(aux_vec,size(boxMEST,1),1),1)/size(boxMEST,1),'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('DEP');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    if(simId==2 && strcmp(sweep_var,'SNR'))
        set(gca,'XTickLabel',sweep_vec-13);
    else
        set(gca,'XTickLabel',sweep_vec);
    end
    grid on;
    %legend(leyenda(1));
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn21/' subcarpeta 'DEP_' sweep_var '_s.pdf']);
    end
    
    % Plot Boxplot
    figure;
    boxplot(boxMEST,'whisker',inf);
    hold on;
    plot(1:length(sweep_vec),mean(boxMEST,1),'o','Color',[1 0.4 1],'MarkerFaceColor',[1 0.4 1])
    ylabel('M_+');
    xlabel(etiquetaX);
    set(gca,'Xtick',1:length(sweep_vec));
    if(simId==2 && strcmp(sweep_var,'SNR'))
        set(gca,'XTickLabel',sweep_vec-13);
    else
        set(gca,'XTickLabel',sweep_vec);
    end
    hold on;
    if(strcmp(sweep_var,'Nt'))
        plot(1:length(sweep_vec),sweep_vec,'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    else
        plot(1:length(sweep_vec),Ntaux*ones(size(sweep_vec)),'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    end
    grid on;
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn21/' subcarpeta 'BoxM_' sweep_var '_s.pdf']);
    end
    
    % Plot Boxplot with #Recovered Cases
    figure;
    boxplot(squeeze(ALL_RECOVERED),'whisker',inf);
    hold on;
    plot(1:length(sweep_vec),squeeze(mean(ALL_RECOVERED,1)),'o','Color',[1 0.4 1],'MarkerFaceColor',[1 0.4 1])
    ylabel('Recovered Tx');
    xlabel(etiquetaX);
    set(gca,'Xtick',1:length(sweep_vec));
    if(simId==2 && strcmp(sweep_var,'SNR'))
        set(gca,'XTickLabel',sweep_vec-13);
    else
        set(gca,'XTickLabel',sweep_vec);
    end
    hold on;
    if(strcmp(sweep_var,'Nt'))
        plot(1:length(sweep_vec),sweep_vec,'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    else
        plot(1:length(sweep_vec),Ntaux*ones(size(sweep_vec)),'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    end
    grid on;
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn21/' subcarpeta 'BoxRecov_' sweep_var '_s.pdf']);
    end
end
