clear all;
close all;
clc

simId = 2;
maxItCluster = 50;
Niter =  10000;

% Base configuration:
T = 500;
Nt = 5;
SNR = 0;
Nr = 20;
L = 1;
Ltrue = 1;
M = 2;

% Sweep and hold variables:
  %sweep_var = 'SNR'; sweep_vec = -5:3:13;
  sweep_var = 'Nt'; sweep_vec = 2:2:10;
  %sweep_var = 'Nr'; sweep_vec = 5:5:25;
  %sweep_var = 'M'; sweep_vec = 2:1:7;
hold_var = 'SNR';
hold_vec = [4];
flagLoadPGAS = 0;
flagLoadBCJR = 0;

% Metrics of interest
ALL_ADER = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MMSE = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_LLH = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));
ALL_MEST = zeros(maxItCluster,Niter+1,length(sweep_vec),length(hold_vec));

ALL_ADER_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_PGAS = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_ADER_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ALL_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));
ALL_SER_ACT_BCJR = zeros(maxItCluster,1,length(sweep_vec),length(hold_vec));

% Load results
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
        else
            error('Not valid value for sweep_var')
        end
        
        % Load all itCluster's for the given configuration
        for itCluster=1:maxItCluster
            
            
            saveFile = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId)...
                        '/T' num2str(Taux) '_Nt' num2str(Ntaux) '_Nr' num2str(Nraux) '_M' num2str(Maux) '_Ltrue' num2str(Ltrueaux) '_L' num2str(Laux) '_SNR' num2str(SNRaux) '/itCluster' num2str(itCluster) '.mat'];
                    
            if(exist(saveFile,'file'))
                load(saveFile,'ADER','MMSE','SER_ACT','SER_ALL','LLH','M_EST');
                ALL_ADER(itCluster,:,idxSweep,idxHold) = ADER;
                ALL_MMSE(itCluster,:,idxSweep,idxHold) = MMSE;
                ALL_LLH(itCluster,:,idxSweep,idxHold) = LLH;
                ALL_SER_ALL(itCluster,:,idxSweep,idxHold) = SER_ALL;
                ALL_SER_ACT(itCluster,:,idxSweep,idxHold) = SER_ACT;
                ALL_MEST(itCluster,:,idxSweep,idxHold) = M_EST;
                if(flagLoadPGAS)
                    load(saveFile,'*_PGAS');
                    if(~exist('ADER_PGAS','var'));
                        idxItClusterNotFound = unique([idxItClusterNotFound itCluster]);
                        disp(['File ' saveFile ' does not contain Genie variables']);
                    end
                elseif(flagLoadBCJR)
                    load(saveFile,'*_BCJR');
                    if(~exist('ADER_BCJR','var'));
                        idxItClusterNotFound = unique([idxItClusterNotFound itCluster]);
                        disp(['File ' saveFile ' does not contain Genie variables']);
                    end
                end
            else
                idxItClusterNotFound = unique([idxItClusterNotFound itCluster]);
                disp(['File ' saveFile ' not found']);
            end
        end
        
        
    end
end

% Remove not found itCluster's
if(~isempty(idxItClusterNotFound))
    disp(['WARNING: Removing ' num2str(length(idxItClusterNotFound)) ' simulations. Total is now: ' num2str(maxItCluster-length(idxItClusterNotFound))]);
end
ALL_ADER(idxItClusterNotFound,:,:,:) = [];
ALL_MMSE(idxItClusterNotFound,:,:,:) = [];
ALL_LLH(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT(idxItClusterNotFound,:,:,:) = [];
ALL_MEST(idxItClusterNotFound,:,:,:) = [];


% For fixed values of the hold variable, plot the metrics vs the sweep variable
for holdV=hold_vec
    idxHold = find(hold_vec==holdV);
    
    figure(idxHold);
    if(strcmp(hold_var,'Nt'))
        Nt = holdV;
    end
    
    % Plot the ADER
    avg_good = zeros(1,length(sweep_vec));
    avg_bad = zeros(1,length(sweep_vec));
    howMany_good = zeros(1,length(sweep_vec));
    howMany_good_text = cell(1,length(sweep_vec));
    subplot(2,3,1);
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Nt = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Nt);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Nt);

        avg_good(idxSweep) = mean(ALL_ADER(idxGood,end,idxSweep,idxHold));
        avg_bad(idxSweep) = mean(ALL_ADER(:,end,idxSweep,idxHold));
        howMany_good(idxSweep) = length(idxGood);
        howMany_good_text{idxSweep} = num2str(length(idxGood));

        plot(sweepV*ones(length(idxGood),1),ALL_ADER(idxGood,end,idxSweep,idxHold),'x','Color','g');
        hold on;
        plot(sweepV*ones(length(idxBad),1),ALL_ADER(idxBad,end,idxSweep,idxHold),'x','Color','r');
    end
    plot(sweep_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(sweep_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',sweep_vec);
    h = text(sweep_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title([hold_var '=' num2str(holdV)]);
    xlabel(sweep_var);
    ylabel('ADER');
    
    % Plot the SER_ALL
    avg_good = zeros(1,length(sweep_vec));
    avg_bad = zeros(1,length(sweep_vec));
    howMany_good = zeros(1,length(sweep_vec));
    howMany_good_text = cell(1,length(sweep_vec));
    subplot(2,3,2);
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Nt = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Nt);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Nt);

        avg_good(idxSweep) = mean(ALL_SER_ALL(idxGood,end,idxSweep,idxHold));
        avg_bad(idxSweep) = mean(ALL_SER_ALL(:,end,idxSweep,idxHold));
        howMany_good(idxSweep) = length(idxGood);
        howMany_good_text{idxSweep} = num2str(length(idxGood));

        plot(sweepV*ones(length(idxGood),1),ALL_SER_ALL(idxGood,end,idxSweep,idxHold),'x','Color','g');
        hold on;
        plot(sweepV*ones(length(idxBad),1),ALL_SER_ALL(idxBad,end,idxSweep,idxHold),'x','Color','r');
    end
    plot(sweep_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(sweep_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',sweep_vec);
    h = text(sweep_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title([hold_var '=' num2str(holdV)]);
    xlabel(sweep_var);
    ylabel('SER(ALL)');

    % Plot the SER_ACT
    avg_good = zeros(1,length(sweep_vec));
    avg_bad = zeros(1,length(sweep_vec));
    howMany_good = zeros(1,length(sweep_vec));
    howMany_good_text = cell(1,length(sweep_vec));
    subplot(2,3,3);
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Nt = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Nt);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Nt);

        avg_good(idxSweep) = mean(ALL_SER_ACT(idxGood,end,idxSweep,idxHold));
        avg_bad(idxSweep) = mean(ALL_SER_ACT(:,end,idxSweep,idxHold));
        howMany_good(idxSweep) = length(idxGood);
        howMany_good_text{idxSweep} = num2str(length(idxGood));

        plot(sweepV*ones(length(idxGood),1),ALL_SER_ACT(idxGood,end,idxSweep,idxHold),'x','Color','g');
        hold on;
        plot(sweepV*ones(length(idxBad),1),ALL_SER_ACT(idxBad,end,idxSweep,idxHold),'x','Color','r');
    end
    plot(sweep_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(sweep_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',sweep_vec);
    h = text(sweep_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title([hold_var '=' num2str(holdV)]);
    xlabel(sweep_var);
    ylabel('SER(ACT)');

    % Plot the LLH
    avg_good = zeros(1,length(sweep_vec));
    avg_bad = zeros(1,length(sweep_vec));
    howMany_good = zeros(1,length(sweep_vec));
    howMany_good_text = cell(1,length(sweep_vec));
    subplot(2,3,4);
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Nt = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Nt);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Nt);

        avg_good(idxSweep) = mean(ALL_LLH(idxGood,end,idxSweep,idxHold));
        avg_bad(idxSweep) = mean(ALL_LLH(:,end,idxSweep,idxHold));
        howMany_good(idxSweep) = length(idxGood);
        howMany_good_text{idxSweep} = num2str(length(idxGood));

        plot(sweepV*ones(length(idxGood),1),ALL_LLH(idxGood,end,idxSweep,idxHold),'x','Color','g');
        hold on;
        plot(sweepV*ones(length(idxBad),1),ALL_LLH(idxBad,end,idxSweep,idxHold),'x','Color','r');
    end
    plot(sweep_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(sweep_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',sweep_vec);
    h = text(sweep_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title([hold_var '=' num2str(holdV)]);
    xlabel(sweep_var);
    ylabel('LLH');
    
    
    % Plot MSE
    avg_good = zeros(1,length(sweep_vec));
    avg_bad = zeros(1,length(sweep_vec));
    howMany_good = zeros(1,length(sweep_vec));
    howMany_good_text = cell(1,length(sweep_vec));
    subplot(2,3,5);
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Nt = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Nt);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Nt);

        avg_good(idxSweep) = mean(ALL_MMSE(idxGood,end,idxSweep,idxHold));
        avg_bad(idxSweep) = mean(ALL_MMSE(:,end,idxSweep,idxHold));
        howMany_good(idxSweep) = length(idxGood);
        howMany_good_text{idxSweep} = num2str(length(idxGood));

        plot(sweepV*ones(length(idxGood),1),ALL_MMSE(idxGood,end,idxSweep,idxHold),'x','Color','g');
        hold on;
        plot(sweepV*ones(length(idxBad),1),ALL_MMSE(idxBad,end,idxSweep,idxHold),'x','Color','r');
    end
    plot(sweep_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(sweep_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',sweep_vec);
    h = text(sweep_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title([hold_var '=' num2str(holdV)]);
    xlabel(sweep_var);
    ylabel('MSE');
    
    % Plot Mest
    avg_good = zeros(1,length(sweep_vec));
    avg_bad = zeros(1,length(sweep_vec));
    howMany_good = zeros(1,length(sweep_vec));
    howMany_good_text = cell(1,length(sweep_vec));
    subplot(2,3,6);
    for sweepV=sweep_vec
        idxSweep = find(sweep_vec==sweepV);
        
        if(strcmp(sweep_var,'Nt'))
            Nt = sweepV;
        end
    
        idxGood = find(ALL_MEST(:,end,idxSweep,idxHold)==Nt);
        idxBad = find(ALL_MEST(:,end,idxSweep,idxHold)~=Nt);

        avg_good(idxSweep) = mean(ALL_MEST(idxGood,end,idxSweep,idxHold));
        avg_bad(idxSweep) = mean(ALL_MEST(:,end,idxSweep,idxHold));
        howMany_good(idxSweep) = length(idxGood);
        howMany_good_text{idxSweep} = num2str(length(idxGood));

        plot(sweepV*ones(length(idxGood),1),ALL_MEST(idxGood,end,idxSweep,idxHold),'x','Color','g');
        hold on;
        plot(sweepV*ones(length(idxBad),1),ALL_MEST(idxBad,end,idxSweep,idxHold),'x','Color','r');
    end
    plot(sweep_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(sweep_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',sweep_vec);
    h = text(sweep_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title([hold_var '=' num2str(holdV)]);
    xlabel(sweep_var);
    ylabel('M_+');
end
