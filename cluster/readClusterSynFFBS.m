clear all;
close all;
clc
addpath aboxplot;

maxItCluster = 50;
Niter =  10000;

%% Scenario definition

%%%%%%%%%%%%%%%
plotToFile = 0;
%%%%%%%%%%%%%%%

%Base scenario
T = 500;
Nt = 5;
SNR = 0;
Nr = 10;
L = 1;
Ltrue = 1;
M = 2;
lHead = 0;
onOffModel = 0;

% Sweep and hold variables:
sweep_var = 'L'; sweep_vec = 1:2; simId = 9; etiquetaX='L'; lugar='NorthEast'; hold_var='SNR'; hold_vec=0;

marcadores = {'+','^','o'};
estilos = {'-','--',':'};
colores = {'m','b','k'};
leyenda = {'PGAS','FFBS'};
anchoSm = 3; altoSm = 2;
%anchoLar = 5; altoLar = 3.5;

%% Initialize metrics of interest
ALL_ADER = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MMSE = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_LLH = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MEST = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));

ALL_ADER_FFBS = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_FFBS = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_FFBS = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MMSE_FFBS = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_LLH_FFBS = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MEST_FFBS = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));

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
    Ltrueaux = L;
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
        Ltrueaux = sweepV;
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
                load(saveFile,'ADER','MMSE','SER_ACT','SER_ALL','LLH','M_EST');
                ALL_ADER(itCluster,:,idxSweep,idxHold) = ADER;
                ALL_MMSE(itCluster,:,idxSweep,idxHold) = MMSE;
                ALL_LLH(itCluster,:,idxSweep,idxHold) = LLH;
                ALL_SER_ALL(itCluster,:,idxSweep,idxHold) = SER_ALL;
                ALL_SER_ACT(itCluster,:,idxSweep,idxHold) = SER_ACT;
                ALL_MEST(itCluster,:,idxSweep,idxHold) = M_EST;
                
                if(exist([saveFile(1:end-4) 'FFBS.mat'],'file'))
                    load([saveFile(1:end-4) 'FFBS.mat'],'ADER','MMSE','SER_ACT','SER_ALL','LLH','M_EST');
                    ALL_ADER_FFBS(itCluster,:,idxSweep,idxHold) = ADER;
                    ALL_SER_ALL_FFBS(itCluster,:,idxSweep,idxHold) = SER_ALL;
                    ALL_SER_ACT_FFBS(itCluster,:,idxSweep,idxHold) = SER_ACT;
                    ALL_MMSE_FFBS(itCluster,:,idxSweep,idxHold) = MMSE;
                    ALL_LLH_FFBS(itCluster,:,idxSweep,idxHold) = LLH;
                    ALL_MEST_FFBS(itCluster,:,idxSweep,idxHold) = M_EST;
                else
                    idxItClusterNotFound = unique([idxItClusterNotFound itCluster]);
                    disp(['File ' saveFile(1:end-4) 'FFBS.mat' ' not found']);
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

ALL_ADER_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_LLH_FFBS(idxItClusterNotFound,:,:,:) = [];
ALL_MEST_FFBS(idxItClusterNotFound,:,:,:) = [];

%% For fixed values of the hold variable, plot the metrics vs the sweep variable
for holdV=hold_vec
    idxHold = find(hold_vec==holdV);
    
    % Go back to base configuration
    Taux = T;
    Ntaux = Nt;
    SNRaux = SNR;
    Nraux = Nr;
    Laux = L;
    Ltrueaux = L;
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
    
    avgADER_FFBS = zeros(1,length(sweep_vec));
    avgSER_ACT_FFBS = zeros(1,length(sweep_vec));
    avgSER_ALL_FFBS = zeros(1,length(sweep_vec));
    avgMSE_FFBS = zeros(1,length(sweep_vec));
    avgMEST_FFBS = zeros(1,length(sweep_vec));
    boxMEST_FFBS = zeros(maxItCluster-length(idxItClusterNotFound),length(sweep_vec));
    
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Ntaux = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Ntaux);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Ntaux);

        idxGood_FFBS = find(ALL_MEST_FFBS(:,end,idxSweep,idxHold)==Ntaux);
        idxBad_FFBS = find(ALL_MEST_FFBS(:,end,idxSweep,idxHold)~=Ntaux);
        
        avgADER(idxSweep) = mean(ALL_ADER(idxGood,end,idxSweep,idxHold));
        avgSER_ACT(idxSweep) = mean(ALL_SER_ACT(idxGood,end,idxSweep,idxHold));
        avgSER_ALL(idxSweep) = mean(ALL_SER_ALL(idxGood,end,idxSweep,idxHold));
        avgMEST(idxSweep) = mean(ALL_MEST(idxGood,end,idxSweep,idxHold));
        avgMSE(idxSweep) = mean(ALL_MMSE(idxGood,end,idxSweep,idxHold));
        boxMEST(:,idxSweep) = ALL_MEST(:,end,idxSweep,idxHold);
        
        avgADER_FFBS(idxSweep) = mean(ALL_ADER_FFBS(idxGood_FFBS,end,idxSweep,idxHold));
        avgSER_ACT_FFBS(idxSweep) = mean(ALL_SER_ACT_FFBS(idxGood_FFBS,end,idxSweep,idxHold));
        avgSER_ALL_FFBS(idxSweep) = mean(ALL_SER_ALL_FFBS(idxGood_FFBS,end,idxSweep,idxHold));
        avgMEST_FFBS(idxSweep) = mean(ALL_MEST_FFBS(idxGood_FFBS,end,idxSweep,idxHold));
        avgMSE_FFBS(idxSweep) = mean(ALL_MMSE_FFBS(idxGood_FFBS,end,idxSweep,idxHold));
        boxMEST_FFBS(:,idxSweep) = ALL_MEST_FFBS(:,end,idxSweep,idxHold);
    end
    
    % Plot ADER
    figure;
    h = plot(sweep_vec,avgADER,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('ADER');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    set(gca,'XTickLabel',sweep_vec);
    grid on;
    hold on;
    h = plot(sweep_vec,avgADER_FFBS,'Marker',marcadores{2},'LineStyle',estilos{2},'Color',colores{2});
    legend(leyenda(1:2),'Location',lugar);
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/ffbs_pgas/ADER_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/ffbs_pgas/ADER_' sweep_var '_g.pdf']);
    end
    
    % Plot SER_ALL
    figure;
    h = plot(sweep_vec,avgSER_ALL,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('SER');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    set(gca,'XTickLabel',sweep_vec);
    hold on;
    grid on;
    h = plot(sweep_vec,avgSER_ALL_FFBS,'Marker',marcadores{2},'LineStyle',estilos{2},'Color',colores{2});
    legend(leyenda(1:2),'Location',lugar);
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/ffbs_pgas/SER_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/ffbs_pgas/SER_' sweep_var '_g.pdf']);
    end
    
    % Plot MSE
    figure;
    h = plot(sweep_vec,avgMSE,'Marker',marcadores{1},'LineStyle',estilos{1},'Color',colores{1});
    ylabel('MSE');
    xlabel(etiquetaX);
    set(gca,'XLim',[min(sweep_vec) max(sweep_vec)]);
    set(gca,'XTick',sweep_vec);
    set(gca,'XTickLabel',sweep_vec);
    hold on;
    grid on;
    h = plot(sweep_vec,avgMSE_FFBS,'Marker',marcadores{2},'LineStyle',estilos{2},'Color',colores{2});
    legend(leyenda(1:2),'Location',lugar);
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/ffbs_pgas/MSE_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/ffbs_pgas/MSE_' sweep_var '_g.pdf']);
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
    set(gca,'XTickLabel',sweep_vec);
    hold on;
    grid on;
    h = plot(sweep_vec,1-sum(boxMEST_FFBS==repmat(aux_vec,size(boxMEST_FFBS,1),1),1)/size(boxMEST_FFBS,1),'Marker',marcadores{2},'LineStyle',estilos{2},'Color',colores{2});
    legend(leyenda(1:2),'Location',lugar);
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/ffbs_pgas/DEP_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/ffbs_pgas/DEP_' sweep_var '_g.pdf']);
    end
    
    % Plot Boxplot
    figure;
    boxplot(boxMEST,'whisker',inf);
    hold on;
    plot(1:length(sweep_vec),mean(boxMEST,1),'o','Color',[1 0.4 1],'MarkerFaceColor',[1 0.4 1])
    ylabel('M_+');
    xlabel(etiquetaX);
    set(gca,'Xtick',1:length(sweep_vec));
    set(gca,'XTickLabel',sweep_vec);
    hold on;
    title(leyenda(1));
    if(strcmp(sweep_var,'Nt'))
        plot(1:length(sweep_vec),sweep_vec,'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    else
        plot(1:length(sweep_vec),Ntaux*ones(size(sweep_vec)),'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    end
    grid on;
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/ffbs_pgas/BoxM_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/ffbs_pgas/BoxM_' sweep_var '_g.pdf']);
    end
    
    % Plot Boxplot (FFBS)
    figure;
    boxplot(boxMEST_FFBS,'whisker',inf);
    hold on;
    plot(1:length(sweep_vec),mean(boxMEST_FFBS,1),'o','Color',[1 0.4 1],'MarkerFaceColor',[1 0.4 1])
    ylabel('M_+');
    xlabel(etiquetaX);
    set(gca,'Xtick',1:length(sweep_vec));
    set(gca,'XTickLabel',sweep_vec);
    hold on;
    title(leyenda(2));
    if(strcmp(sweep_var,'Nt'))
        plot(1:length(sweep_vec),sweep_vec,'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    else
        plot(1:length(sweep_vec),Ntaux*ones(size(sweep_vec)),'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
    end
    grid on;
    if(plotToFile)
        figurapdf(anchoSm,altoSm);
        print('-dpdf',['./plots/ffbs_pgas/BoxM_' sweep_var '_s.pdf']);
        figurapdf(anchoLar,altoLar);
        print('-dpdf',['./plots/ffbs_pgas/BoxM_' sweep_var '_g.pdf']);
    end    
    
end


% Plot all boxplots
figure;
anchoLar = 4.5; altoLar = 3;
boxAUX = cat(1,permute(boxMEST_FFBS,[3 1 2]),permute(boxMEST,[3 1 2]));
aboxplot(boxAUX,'labels',sweep_vec);
set(gca,'FontSize',14);
legend('FFBS','PGAS','Location','NorthWest');
xlabel(etiquetaX);
ylabel('M_+');
hold on;
set(gca,'yLim',[min(boxAUX(:))-0.5,max(boxAUX(:))+0.5]);
set(gca,'YGrid','on');
set(gca,'YTick',[min(boxAUX(:)):1:max(boxAUX(:))]);
deltaXaxis = 0.175;
hold on;
if(strcmp(sweep_var,'Nt'))
    plot([1-deltaXaxis 1+deltaXaxis 2-deltaXaxis 2+deltaXaxis],sweep_vec,'p','Color',[0 0.4 0],'MarkerFaceColor',[0 0.4 0],'MarkerSize',10);
else
    plot([1-deltaXaxis 1+deltaXaxis 2-deltaXaxis 2+deltaXaxis],Ntaux*ones(1,4),'p','Color',[0 0.4 0],'MarkerFaceColor',[0 0.4 0],'MarkerSize',10);
end
if(plotToFile)
    figurapdf(anchoSm,altoSm);
    print('-dpdf',['./plots/ffbs_pgas/BoxM_' sweep_var '_s.pdf']);
    figurapdf(anchoLar,altoLar);
    print('-dpdf',['./plots/ffbs_pgas/BoxM_' sweep_var '_g.pdf']);
end    

