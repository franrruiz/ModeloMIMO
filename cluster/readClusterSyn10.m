clear all;
close all;
clc
addpath aboxplot;

maxItCluster = 50;
Niter =  10000;

%% Scenario definition


%%%%%%%%%%%%%%%
pRemove = 0.0;
plotToFile = 1;
flagPlotGenie = 0;
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

% Sweep and hold variables:
  % L = 1:
  %sweep_var = 'lHead'; sweep_vec = 0:3:3; simId=10; etiquetaX='Header Length'; lugar='NorthEast'; hold_var='onOffModel'; hold_vec=0:1;
  %sweep_var = 'SNR'; sweep_vec = -18:3:0; simId = 10; etiquetaX='SNR (dB)'; lugar='NorthEast'; hold_var='Nr'; hold_vec=20;
  %sweep_var = 'Nt'; sweep_vec = 2:2:10; simId = 10; etiquetaX='#Transmitters'; lugar='NorthWest'; hold_var='SNR'; hold_vec=-3;
  %sweep_var = 'Nr'; sweep_vec = 5:5:20; simId = 10; etiquetaX='#Receivers'; lugar='NorthEast'; hold_var='SNR'; hold_vec=-3;
  %sweep_var = 'M'; sweep_vec = 2:1:7; simId = 10; etiquetaX='log_2|A|'; lugar='NorthWest'; hold_var='SNR'; hold_vec=-3;

  % L = 3:
  sweep_var = 'SNR'; sweep_vec = [-18 -15 -12 -3:3:0]; simId = 10; etiquetaX='SNR (dB)'; lugar='NorthEast'; hold_var='L'; hold_vec=3; Ltrue=3;
  %sweep_var = 'M'; sweep_vec = 2:2:4; simId = 10; etiquetaX='log_2|A|'; lugar='NorthWest'; hold_var='L'; hold_vec=3; Ltrue=3;
  %sweep_var = 'SNR'; sweep_vec = -18:3:-12; simId = 10; etiquetaX='SNR (dB)'; lugar='NorthEast'; hold_var='L'; hold_vec=2;

marcadores = {'+','^','o'};
estilos = {'-','--',':'};
colores = {'m','b','k'};
leyenda = {'iFHMM','PGAS','BCJR'};
anchoSm = 3; altoSm = 2;
anchoLar = 5; altoLar = 3.5;

%% Initialize metrics of interest
ALL_ADER = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MMSE = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_LLH = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MEST = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));

ALL_ADER_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_MMSE_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_ADER_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_MMSE_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));

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
                        '/T' num2str(Taux) '_Nt' num2str(Ntaux) '_Nr' num2str(Nraux) '_M' num2str(Maux) '_Ltrue' num2str(Ltrueaux) '_L' num2str(Laux) '_SNR' num2str(SNRaux) '_lHead' num2str(lHeadaux), '_onOff' num2str(onOffModelaux) ...
                        '/itCluster' num2str(itCluster) '.mat'];
                    
            if(exist(saveFile,'file'))
                clear ADER_PGAS2 ADER_BCJR2
                load(saveFile,'ADER','MMSE','SER_ACT','SER_ALL','LLH','M_EST','*_PGAS2','*_BCJR2');
                ALL_ADER(itCluster,:,idxSweep,idxHold) = ADER;
                ALL_MMSE(itCluster,:,idxSweep,idxHold) = MMSE;
                ALL_LLH(itCluster,:,idxSweep,idxHold) = LLH;
                ALL_SER_ALL(itCluster,:,idxSweep,idxHold) = SER_ALL;
                ALL_SER_ACT(itCluster,:,idxSweep,idxHold) = SER_ACT;
                ALL_MEST(itCluster,:,idxSweep,idxHold) = M_EST;
                
                if(exist('ADER_PGAS2','var'))
                    ALL_ADER_PGAS(itCluster,1,idxSweep,idxHold) = ADER_PGAS2;
                    ALL_SER_ALL_PGAS(itCluster,1,idxSweep,idxHold) = SER_ALL_PGAS2;
                    ALL_SER_ACT_PGAS(itCluster,1,idxSweep,idxHold) = SER_ACT_PGAS2;
                    ALL_MMSE_PGAS(itCluster,1,idxSweep,idxHold) = MMSE_PGAS2;
                end
                if(exist('ADER_BCJR2','var'))
                    ALL_ADER_BCJR(itCluster,1,idxSweep,idxHold) = ADER_BCJR2;
                    ALL_SER_ALL_BCJR(itCluster,1,idxSweep,idxHold) = SER_ALL_BCJR2;
                    ALL_SER_ACT_BCJR(itCluster,1,idxSweep,idxHold) = SER_ACT_BCJR2;
                    ALL_MMSE_BCJR(itCluster,1,idxSweep,idxHold) = MMSE_BCJR2;
                end
            else
                idxItClusterNotFound = unique([idxItClusterNotFound itCluster]);
                disp(['File ' saveFile ' not found']);
            end
        end
        
    end
end

%% Remove not found itCluster's
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


%% Remove pRemove% of the worst simulations
Ntaux = Nt;
for idxSweep=1:length(sweep_vec)
    if(strcmp(sweep_var,'Nt'))
        Ntaux = sweep_vec(idxSweep);
    end
    for idxHold=1:length(hold_vec)
        if(strcmp(hold_var,'Nt'))
            Ntaux = hold_vec(idxHold);
        end
        
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Ntaux);
        howManyToRm = round(pRemove*length(idxGood));
        [valnul idxOrd] = sort(ALL_MMSE(idxGood,end,idxSweep,idxHold),'descend');
        ALL_MEST(idxGood(idxOrd(1:howManyToRm)),end,idxSweep,idxHold) = -ALL_MEST(idxGood(idxOrd(1:howManyToRm)),end,idxSweep,idxHold);
    end
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
    
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Ntaux = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Ntaux);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Ntaux);
        
        avgADER(idxSweep) = mean(ALL_ADER(idxGood,end,idxSweep,idxHold));
        stdADER(idxSweep) = std(ALL_ADER(idxGood,end,idxSweep,idxHold));
        avgSER_ACT(idxSweep) = mean(ALL_SER_ACT(idxGood,end,idxSweep,idxHold));
        stdSER_ACT(idxSweep) = std(ALL_SER_ACT(idxGood,end,idxSweep,idxHold));
        avgSER_ALL(idxSweep) = mean(ALL_SER_ALL(idxGood,end,idxSweep,idxHold));
        stdSER_ALL(idxSweep) = std(ALL_SER_ALL(idxGood,end,idxSweep,idxHold));
        avgMEST(idxSweep) = mean(ALL_MEST(idxGood,end,idxSweep,idxHold));
        stdMEST(idxSweep) = std(ALL_MEST(idxGood,end,idxSweep,idxHold));
        avgMSE(idxSweep) = mean(ALL_MMSE(idxGood,end,idxSweep,idxHold));
        stdMSE(idxSweep) = std(ALL_MMSE(idxGood,end,idxSweep,idxHold));
        boxMEST(:,idxSweep) = abs(ALL_MEST(:,end,idxSweep,idxHold));
        
        avgADER_PGAS(idxSweep) = mean(ALL_ADER_PGAS(idxGood,1,idxSweep,idxHold));
        stdADER_PGAS(idxSweep) = std(ALL_ADER_PGAS(idxGood,1,idxSweep,idxHold));
        avgSER_ACT_PGAS(idxSweep) = mean(ALL_SER_ACT_PGAS(idxGood,1,idxSweep,idxHold));
        stdSER_ACT_PGAS(idxSweep) = std(ALL_SER_ACT_PGAS(idxGood,1,idxSweep,idxHold));
        avgSER_ALL_PGAS(idxSweep) = mean(ALL_SER_ALL_PGAS(idxGood,1,idxSweep,idxHold));
        stdSER_ALL_PGAS(idxSweep) = std(ALL_SER_ALL_PGAS(idxGood,1,idxSweep,idxHold));
        
        avgADER_BCJR(idxSweep) = mean(ALL_ADER_BCJR(idxGood,1,idxSweep,idxHold));
        stdADER_BCJR(idxSweep) = std(ALL_ADER_BCJR(idxGood,1,idxSweep,idxHold));
        avgSER_ACT_BCJR(idxSweep) = mean(ALL_SER_ACT_BCJR(idxGood,1,idxSweep,idxHold));
        stdSER_ACT_BCJR(idxSweep) = std(ALL_SER_ACT_BCJR(idxGood,1,idxSweep,idxHold));
        avgSER_ALL_BCJR(idxSweep) = mean(ALL_SER_ALL_BCJR(idxGood,1,idxSweep,idxHold));
        stdSER_ALL_BCJR(idxSweep) = std(ALL_SER_ALL_BCJR(idxGood,1,idxSweep,idxHold));
    end
    
    % Plot ADER
    figure;
    h = plot(sweep_vec,avgADER,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('ADER');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
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
    if(strcmp(sweep_var,'Nt'))
        idxBCJR = find(((1+2^Maux).^(2*Laux*sweep_vec)<1e6)&(sweep_vec<=4));
    elseif(strcmp(sweep_var,'M'))
        idxBCJR = find(((1+2.^sweep_vec).^(2*Laux*Ntaux)<1e6)&(Ntaux<=4));
    elseif(strcmp(sweep_var,'L'))
        idxBCJR = find(((1+2^Maux).^(2*sweep_vec*Ntaux)<1e6)&(Ntaux<=4));
    else
        idxBCJR = find(((1+2^Maux)^(2*Laux*Ntaux)<1e6)&(Ntaux<=4));
    end
    if(~isempty(idxBCJR))
        hold on;
        h = plot(sweep_vec(idxBCJR),avgADER_BCJR(idxBCJR),'Marker',marcadores{3},'LineStyle',estilos{3},'Color',colores{3});
        if(~flagPlotGenie)
            set(h,'Visible','off');
        else
            legend(leyenda,'Location',lugar);
        end
    elseif(flagPlotGenie)
        legend(leyenda(1:2),'Location',lugar);
    end
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn10/ADER_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/syn10/ADER_' sweep_var '_g.pdf']);
    end
    
    % Plot SER_ALL
    figure;
    h = plot(sweep_vec,avgSER_ALL,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('SER');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
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
    idxBCJR = [];
    if(strcmp(sweep_var,'Nt'))
        idxBCJR = find(((1+2^Maux).^(2*Laux*sweep_vec)<1e6)&(sweep_vec<=4));
    elseif(strcmp(sweep_var,'M'))
        idxBCJR = find(((1+2.^sweep_vec).^(2*Laux*Ntaux)<1e6)&(Ntaux<=4));
    elseif(strcmp(sweep_var,'L'))
        idxBCJR = find(((1+2^Maux).^(2*sweep_vec*Ntaux)<1e6)&(Ntaux<=4));
    else
        idxBCJR = find(((1+2^Maux)^(2*Laux*Ntaux)<1e6)&(Ntaux<=4));
    end
    if(~isempty(idxBCJR))
        hold on;
        h = plot(sweep_vec(idxBCJR),avgSER_ALL_BCJR(idxBCJR),'Marker',marcadores{3},'LineStyle',estilos{3},'Color',colores{3});
        if(~flagPlotGenie)
            set(h,'Visible','off');
        else
            legend(leyenda,'Location',lugar);
        end
    elseif(flagPlotGenie)
        legend(leyenda(1:2),'Location',lugar);
    end
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/syn10/SER_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/syn10/SER_' sweep_var '_g.pdf']);
    end
    
    % Plot MSE
    figure;
    h = plot(sweep_vec,avgMSE,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('MSE');
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
        print('-dpdf',['./plots/syn10/MSE_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/syn10/MSE_' sweep_var '_g.pdf']);
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
        print('-dpdf',['./plots/syn10/DEP_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/syn10/DEP_' sweep_var '_g.pdf']);
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
        print('-dpdf',['./plots/syn10/BoxM_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/syn10/BoxM_' sweep_var '_g.pdf']);
    end
end
