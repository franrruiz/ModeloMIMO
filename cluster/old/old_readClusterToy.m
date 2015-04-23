clear all;
close all;
clc

toyType = 1;
subToyType = 1;
Nt_vec = 2;%:5;
SNR_vec = -20:5:5;%[-10:10:20 50];
maxItCluster = 50;
Niter =  5000;
saveFolder = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/toy' num2str(toyType) '/sub' num2str(subToyType)];

ALL_ADER = zeros(maxItCluster,Niter,length(Nt_vec),length(SNR_vec));
ALL_LLH = zeros(maxItCluster,Niter,length(Nt_vec),length(SNR_vec));
ALL_SER_ALL = zeros(maxItCluster,Niter,length(Nt_vec),length(SNR_vec));
ALL_SER_ACT = zeros(maxItCluster,Niter,length(Nt_vec),length(SNR_vec));
ALL_MEST = zeros(maxItCluster,Niter,length(Nt_vec),length(SNR_vec));

idxItClusterNotFound = [];
for Nt=Nt_vec
    idxNt = find(Nt==Nt_vec);
    for SNR=SNR_vec
        idxSNR = find(SNR==SNR_vec);
        for itCluster=1:maxItCluster
            saveFile = [saveFolder '/Nt' num2str(Nt) '_SNR' num2str(SNR) '_itCluster' num2str(itCluster) '.mat'];
            
            if(exist(saveFile,'file'))
                load(saveFile,'ADER','SER_ACT','SER_ALL','LLH','M_EST');

                ALL_ADER(itCluster,:,idxNt,idxSNR) = ADER;
                ALL_LLH(itCluster,:,idxNt,idxSNR) = LLH;
                ALL_SER_ALL(itCluster,:,idxNt,idxSNR) = SER_ALL;
                ALL_SER_ACT(itCluster,:,idxNt,idxSNR) = SER_ACT;
                ALL_MEST(itCluster,:,idxNt,idxSNR) = M_EST;
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
ALL_LLH(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ALL(idxItClusterNotFound,:,:,:) = [];
ALL_SER_ACT(idxItClusterNotFound,:,:,:) = [];
ALL_MEST(idxItClusterNotFound,:,:,:) = [];

%% For fixed Nt: Plot ADER, SER_ALL, SER_ACT vs SNR
figure(1);
for Nt=Nt_vec
    idxNt = find(Nt==Nt_vec);
    
    % ADER
    subplot(length(Nt_vec),4,(idxNt-1)*4+1);
    avg_good = zeros(1,length(SNR_vec));
    avg_bad = zeros(1,length(SNR_vec));
    howMany_good = zeros(1,length(SNR_vec));
    howMany_good_text = cell(1,length(SNR_vec));
    for SNR=SNR_vec
        idxSNR = find(SNR==SNR_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxSNR) = mean(ALL_ADER(idxGood,end,idxNt,idxSNR));
        avg_bad(idxSNR) = mean(ALL_ADER(:,end,idxNt,idxSNR));
        howMany_good(idxSNR) = length(idxGood);
        howMany_good_text{idxSNR} = num2str(length(idxGood));

        plot(SNR*ones(length(idxGood),1),ALL_ADER(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(SNR*ones(length(idxBad),1),ALL_ADER(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(SNR_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(SNR_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',SNR_vec);
    h = text(SNR_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['Nt=' num2str(Nt)]);
    xlabel('SNR (dB)');
    ylabel('ADER');
    
    % SER_ALL
    subplot(length(Nt_vec),4,(idxNt-1)*4+2);
    avg_good = zeros(1,length(SNR_vec));
    avg_bad = zeros(1,length(SNR_vec));
    howMany_good = zeros(1,length(SNR_vec));
    howMany_good_text = cell(1,length(SNR_vec));
    for SNR=SNR_vec
        idxSNR = find(SNR==SNR_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxSNR) = mean(ALL_SER_ALL(idxGood,end,idxNt,idxSNR));
        avg_bad(idxSNR) = mean(ALL_SER_ALL(:,end,idxNt,idxSNR));
        howMany_good(idxSNR) = length(idxGood);
        howMany_good_text{idxSNR} = num2str(length(idxGood));

        plot(SNR*ones(length(idxGood),1),ALL_SER_ALL(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(SNR*ones(length(idxBad),1),ALL_SER_ALL(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(SNR_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(SNR_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',SNR_vec);
    h = text(SNR_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['Nt=' num2str(Nt)]);
    xlabel('SNR (dB)');
    ylabel('SER (ALL)');
    
    % SER_ACT
    subplot(length(Nt_vec),4,(idxNt-1)*4+3);
    avg_good = zeros(1,length(SNR_vec));
    avg_bad = zeros(1,length(SNR_vec));
    howMany_good = zeros(1,length(SNR_vec));
    howMany_good_text = cell(1,length(SNR_vec));
    for SNR=SNR_vec
        idxSNR = find(SNR==SNR_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxSNR) = mean(ALL_SER_ACT(idxGood,end,idxNt,idxSNR));
        avg_bad(idxSNR) = mean(ALL_SER_ACT(:,end,idxNt,idxSNR));
        howMany_good(idxSNR) = length(idxGood);
        howMany_good_text{idxSNR} = num2str(length(idxGood));

        plot(SNR*ones(length(idxGood),1),ALL_SER_ACT(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(SNR*ones(length(idxBad),1),ALL_SER_ACT(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(SNR_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(SNR_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',SNR_vec);
    h = text(SNR_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['Nt=' num2str(Nt)]);
    xlabel('SNR (dB)');
    ylabel('SER (ACT)');

    % LLH
    subplot(length(Nt_vec),4,(idxNt-1)*4+4);
    avg_good = zeros(1,length(SNR_vec));
    avg_bad = zeros(1,length(SNR_vec));
    howMany_good = zeros(1,length(SNR_vec));
    howMany_good_text = cell(1,length(SNR_vec));
    for SNR=SNR_vec
        idxSNR = find(SNR==SNR_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxSNR) = mean(ALL_LLH(idxGood,end,idxNt,idxSNR));
        avg_bad(idxSNR) = mean(ALL_LLH(:,end,idxNt,idxSNR));
        howMany_good(idxSNR) = length(idxGood);
        howMany_good_text{idxSNR} = num2str(length(idxGood));

        plot(SNR*ones(length(idxGood),1),ALL_LLH(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(SNR*ones(length(idxBad),1),ALL_LLH(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(SNR_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(SNR_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',SNR_vec);
    h = text(SNR_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['Nt=' num2str(Nt)]);
    xlabel('SNR (dB)');
    ylabel('LLH');     
end

%% For fixed SNR: Plot ADER, SER_ALL, SER_ACT vs Nt
figure(2);
for SNR=SNR_vec
    idxSNR = find(SNR==SNR_vec);
    
    % ADER
    subplot(length(SNR_vec),4,(idxSNR-1)*4+1);
    avg_good = zeros(1,length(Nt_vec));
    avg_bad = zeros(1,length(Nt_vec));
    howMany_good = zeros(1,length(Nt_vec));
    howMany_good_text = cell(1,length(Nt_vec));
    for Nt=Nt_vec
        idxNt = find(Nt==Nt_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxNt) = mean(ALL_ADER(idxGood,end,idxNt,idxSNR));
        avg_bad(idxNt) = mean(ALL_ADER(:,end,idxNt,idxSNR));
        howMany_good(idxNt) = length(idxGood);
        howMany_good_text{idxNt} = num2str(length(idxGood));

        plot(Nt*ones(length(idxGood),1),ALL_ADER(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(Nt*ones(length(idxBad),1),ALL_ADER(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(Nt_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(Nt_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',Nt_vec);
    h = text(Nt_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['SNR = ' num2str(SNR) ' dB']);
    xlabel('Nt');
    ylabel('ADER');
    
    % SER_ALL
    subplot(length(SNR_vec),4,(idxSNR-1)*4+2);
    avg_good = zeros(1,length(Nt_vec));
    avg_bad = zeros(1,length(Nt_vec));
    howMany_good = zeros(1,length(Nt_vec));
    howMany_good_text = cell(1,length(Nt_vec));
    for Nt=Nt_vec
        idxNt = find(Nt==Nt_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxNt) = mean(ALL_SER_ALL(idxGood,end,idxNt,idxSNR));
        avg_bad(idxNt) = mean(ALL_SER_ALL(:,end,idxNt,idxSNR));
        howMany_good(idxNt) = length(idxGood);
        howMany_good_text{idxNt} = num2str(length(idxGood));

        plot(Nt*ones(length(idxGood),1),ALL_SER_ALL(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(Nt*ones(length(idxBad),1),ALL_SER_ALL(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(Nt_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(Nt_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',Nt_vec);
    h = text(Nt_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['SNR = ' num2str(SNR) ' dB']);
    xlabel('Nt');
    ylabel('SER (ALL)');
    
    % SER_ACT
    subplot(length(SNR_vec),4,(idxSNR-1)*4+3);
    avg_good = zeros(1,length(Nt_vec));
    avg_bad = zeros(1,length(Nt_vec));
    howMany_good = zeros(1,length(Nt_vec));
    howMany_good_text = cell(1,length(Nt_vec));
    for Nt=Nt_vec
        idxNt = find(Nt==Nt_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxNt) = mean(ALL_SER_ACT(idxGood,end,idxNt,idxSNR));
        avg_bad(idxNt) = mean(ALL_SER_ACT(:,end,idxNt,idxSNR));
        howMany_good(idxNt) = length(idxGood);
        howMany_good_text{idxNt} = num2str(length(idxGood));

        plot(Nt*ones(length(idxGood),1),ALL_SER_ACT(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(Nt*ones(length(idxBad),1),ALL_SER_ACT(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(Nt_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(Nt_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',Nt_vec);
    h = text(Nt_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['SNR = ' num2str(SNR) ' dB']);
    xlabel('Nt');
    ylabel('SER (ACT)');
    
    % LLH
    subplot(length(SNR_vec),4,(idxSNR-1)*4+4);
    avg_good = zeros(1,length(Nt_vec));
    avg_bad = zeros(1,length(Nt_vec));
    howMany_good = zeros(1,length(Nt_vec));
    howMany_good_text = cell(1,length(Nt_vec));
    for Nt=Nt_vec
        idxNt = find(Nt==Nt_vec);
        idxGood = find(ALL_MEST(:,end,idxNt,idxSNR)==Nt);
        idxBad = find(ALL_MEST(:,end,idxNt,idxSNR)~=Nt);

        avg_good(idxNt) = mean(ALL_LLH(idxGood,end,idxNt,idxSNR));
        avg_bad(idxNt) = mean(ALL_LLH(:,end,idxNt,idxSNR));
        howMany_good(idxNt) = length(idxGood);
        howMany_good_text{idxNt} = num2str(length(idxGood));

        plot(Nt*ones(length(idxGood),1),ALL_LLH(idxGood,end,idxNt,idxSNR),'x','Color','g');
        hold on;
        plot(Nt*ones(length(idxBad),1),ALL_LLH(idxBad,end,idxNt,idxSNR),'x','Color','r');
    end
    plot(Nt_vec,avg_good,'Marker','+','LineWidth',1.5,'Color','g');
    hold on;
    plot(Nt_vec,avg_bad,'Marker','+','LineWidth',1.5,'Color','r');
    set(gca,'XTick',Nt_vec);
    h = text(Nt_vec,avg_good,howMany_good_text);
    set(h,'Color','b');
    title(['SNR = ' num2str(SNR) ' dB']);
    xlabel('Nt');
    ylabel('LLH');       
end


%% For fixed Nt: Plot avg M vs SNR
figure(3);
plot(SNR_vec,squeeze(mean(ALL_MEST(:,end,:,:),1)),'Marker','+');
set(gca,'XTick',SNR_vec);
xlabel('SNR (dB)');
ylabel('Avg M_+');
leyenda = cell(1,length(Nt_vec));
for idxNt=1:length(Nt_vec)
    leyenda{idxNt} = ['Nt=' num2str(Nt_vec(idxNt))];
end
for Nt=Nt_vec
    idxNt = find(Nt==Nt_vec);
    for SNR=SNR_vec
        idxSNR = find(SNR==SNR_vec);
        text(SNR,mean(ALL_MEST(:,end,idxNt,idxSNR)),num2str(sum(ALL_MEST(:,end,idxNt,idxSNR)==Nt)));
    end
end
legend(leyenda,'Location','EastOutside');
title('Average value of M_+ for the different scenarios');
grid on;
