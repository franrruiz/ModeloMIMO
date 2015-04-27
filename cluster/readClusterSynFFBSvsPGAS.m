% FFBS vs PGAS
close all; clear all; clc;
addpath aboxplot;

%%%%%%%%%%%%%%%
plotToFile = 1;
%%%%%%%%%%%%%%%

%% Parameters of the scenario to be loaded
maxItCluster = 50;
T = 500;
Nt = 5;
Nr = 10;
M = 2;
SNR = 0;
lHead = 0;
onOffModel = 0;
Nparticles = 3000;
simId = 22;   % simId=22 for these experiments

maxNiter = 10000;
thrSER = 0.1;
L_vec = 1:3;

%% Vary the number of particles
ALL_BOXREC = zeros(maxItCluster,length(L_vec));
ALL_BOXREC_FFBS = zeros(maxItCluster,length(L_vec));
ALL_MEST = zeros(maxItCluster,length(L_vec));
ALL_MEST_FFBS = zeros(maxItCluster,length(L_vec));
for L = L_vec
    Ltrue = L;
    c = find(L_vec==L);
    for itCluster = 1:maxItCluster
    	saveFile = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId)...
                                '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_Ltrue' num2str(Ltrue) '_L' num2str(L) '_SNR' num2str(SNR) '_lHead' num2str(lHead), '_onOff' num2str(onOffModel) '_Npart' num2str(Nparticles) ...
                                '/itCluster' num2str(itCluster) '.mat'];
                        
        load(saveFile,'SER_ALL_indiv','SER_ALL_FFBS_indiv','M_EST','M_EST_FFBS');
        
        ALL_BOXREC(itCluster,c) = sum(SER_ALL_indiv<thrSER);
        ALL_BOXREC_FFBS(itCluster,c) = sum(SER_ALL_FFBS_indiv<thrSER);

        ALL_MEST(itCluster,c) = M_EST(maxNiter+1);
        ALL_MEST_FFBS(itCluster,c) = M_EST_FFBS(maxNiter+1);
    end
    
end

%% M_EST
figure;
boxAUX = cat(1,permute(ALL_MEST_FFBS,[3 1 2]),permute(ALL_MEST,[3 1 2]));
aboxplot(boxAUX,'labels',L_vec);
hold on;
set(gca,'FontSize',14);
legend('FFBS','PGAS','Location','NorthWest');
xlabel('L');
ylabel('M_+');
hold on;
set(gca,'yLim',[min(boxAUX(:))-0.5,max(boxAUX(:))+0.5]);
set(gca,'YGrid','on');
set(gca,'YTick',[min(boxAUX(:)):1:max(boxAUX(:))]);
deltaXaxis = 0.175;
hold on;
plot([1-deltaXaxis 1+deltaXaxis 2-deltaXaxis 2+deltaXaxis 3-deltaXaxis 3+deltaXaxis 4-deltaXaxis 4+deltaXaxis],Nt*ones(1,8),'p','Color',[0 0.4 0],'MarkerFaceColor',[0 0.4 0],'MarkerSize',10);
if(plotToFile)
    figurapdf(4.5,3);  % Before: 3,2
    print('-dpdf',['./plots/syn21/ffbs_pgas/Mest.pdf']);
end

%% M_EST
figure;
boxAUX = cat(1,permute(ALL_BOXREC_FFBS,[3 1 2]),permute(ALL_BOXREC,[3 1 2]));
aboxplot(boxAUX,'labels',L_vec);
hold on;
set(gca,'FontSize',14);
legend('FFBS','PGAS','Location','SouthWest');
xlabel('L');
ylabel('Recovered Tx');
hold on;
set(gca,'yLim',[min(boxAUX(:))-0.5,max(boxAUX(:))+0.5]);
set(gca,'YGrid','on');
set(gca,'YTick',[min(boxAUX(:)):1:max(boxAUX(:))]);
deltaXaxis = 0.175;
hold on;
plot([1-deltaXaxis 1+deltaXaxis 2-deltaXaxis 2+deltaXaxis 3-deltaXaxis 3+deltaXaxis 4-deltaXaxis 4+deltaXaxis],Nt*ones(1,8),'p','Color',[0 0.4 0],'MarkerFaceColor',[0 0.4 0],'MarkerSize',10);
if(plotToFile)
    figurapdf(4,3);  % Before: 3,2
    print('-dpdf',['./plots/syn21/ffbs_pgas/BoxRecov.pdf']);
end
