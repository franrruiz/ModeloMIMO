close all; clear all; clc;

%%%%%%%%%%%%%%%
plotToFile = 1;
%%%%%%%%%%%%%%%

%% Parameters of the scenario to be loaded
itCluster = 1;
T = 1000;
Nt = 10;
Nr = 20;
M = 2;
Ltrue = 1;
L = 1;
SNR = -3;
lHead = 0;
onOffModel = 0;
simId = 21;   % simId=21 for these experiments

maxNiter = 10000;
thrSER = 0.1;

%% Vary the number of particles
part_vec = [300 1000 3000 10000 30000];
ALL_LLH = zeros(length(part_vec),maxNiter);
ALL_SER_ALL_indiv = cell(1,length(part_vec));
ALL_SER_ALL_PGAS_indiv = cell(1,length(part_vec));
ALL_SER_ALL_FFBS_indiv = cell(1,length(part_vec));
ALL_RECOV = zeros(3,length(part_vec));
ALL_MEST = zeros(1,length(part_vec));
for Nparticles = part_vec
    c = find(part_vec==Nparticles);
	saveFile = ['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/' num2str(simId)...
                            '/T' num2str(T) '_Nt' num2str(Nt) '_Nr' num2str(Nr) '_M' num2str(M) '_Ltrue' num2str(Ltrue) '_L' num2str(L) '_SNR' num2str(SNR) '_lHead' num2str(lHead), '_onOff' num2str(onOffModel) '_Npart' num2str(Nparticles) ...
                            '/itCluster' num2str(itCluster) '.mat'];
                        
    load(saveFile,'LLH','SER_ALL_indiv','SER_ALL_PGAS4_indiv','SER_ALL_FFBS4_indiv','M_EST');
    
    ALL_LLH(c,:) = LLH(1:maxNiter);
    ALL_SER_ALL_indiv{c} = SER_ALL_indiv(SER_ALL_indiv<thrSER);
    ALL_SER_ALL_PGAS_indiv{c} = SER_ALL_PGAS4_indiv(SER_ALL_PGAS4_indiv<thrSER);
    ALL_SER_ALL_FFBS_indiv{c} = SER_ALL_FFBS4_indiv(SER_ALL_FFBS4_indiv<thrSER);
    
    ALL_RECOV(1,c) = length(ALL_SER_ALL_indiv{c});
    ALL_RECOV(2,c) = length(ALL_SER_ALL_PGAS_indiv{c});
    ALL_RECOV(3,c) = length(ALL_SER_ALL_FFBS_indiv{c});
    
    ALL_MEST(c) = (M_EST(end));
end

%% Plot LLH vs iterations
figure;
h = plot(1:maxNiter,ALL_LLH);
leyenda = cell(1,length(part_vec));
for i=1:length(part_vec)
    leyenda{i} = [num2str(part_vec(i)) ' particles'];
end
set(h(1),'LineStyle','-','Color',[0 0.6 0]);
set(h(2),'LineStyle','-.','Color','m');
set(h(3),'LineStyle','-','Color',[1 0.5 0]);
set(h(4),'LineStyle','--','Color','r');
set(h(5),'LineStyle',':','Color','k');
legend(leyenda,'Location','SouthEast');
grid on;
xlabel('Iterations');
ylabel('log-likelihood');
set(gca,'YLim',[-8 -5.5]*1e4);
if(plotToFile)
    figurapdf(4,3);  % Before: 3,2
    print('-dpdf',['./plots/syn21/Npart/LLH_Npart_s.pdf']);
end

%% Plot SER vs Npart
figure;
avgSER = zeros(3,length(part_vec));
for c=1:length(part_vec)
    avgSER(1,c) = mean(ALL_SER_ALL_indiv{c});
    avgSER(2,c) = mean(ALL_SER_ALL_PGAS_indiv{c});
    avgSER(3,c) = mean(ALL_SER_ALL_FFBS_indiv{c});
end
plot(part_vec,avgSER,'Marker','+');
set(gca,'XScale','log');
%set(gca,'XTick',part_vec);
%set(gca,'XTickLabel',part_vec);
xlabel('Number of particles');
ylabel('SER');
grid on;
legend('IFFSM','PGAS','FFBS');
if(plotToFile)
    figurapdf(4,3);  % Before: 3,2
    print('-dpdf',['./plots/syn21/Npart/SER_Npart_s.pdf']);
end

%% Number of recovered transmitters
figure;
bar(1:length(part_vec),ALL_RECOV');
hold on;
plot(1:length(part_vec),Nt*ones(1,length(part_vec)),'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
grid on;
set(gca,'XTickLabel',part_vec);
set(gca,'YLim',[0 11]);
xlabel('Number of particles');
ylabel('Recovered Tx');
legend('IFFSM','PGAS','FFBS','Location','EastOutside');
if(plotToFile)
    figurapdf(7,3);  % Before: 3,2
    print('-dpdf',['./plots/syn21/Npart/BoxRecov_Npart_s.pdf']);
end

%% M_EST
figure;
semilogx(part_vec,ALL_MEST,'-','Marker','^');
hold on;
semilogx(part_vec,Nt*ones(1,length(part_vec)),'p','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8);
grid on;
%set(gca,'XTick',part_vec);
%set(gca,'XTickLabel',part_vec);
set(gca,'YLim',[6 11]);
xlabel('Number of particles');
ylabel('M_+');
%legend('IFFSM','PGAS','FFBS','Location','EastOutside');
if(plotToFile)
    figurapdf(4,3);  % Before: 3,2
    print('-dpdf',['./plots/syn21/Npart/Mest_Npart_s.pdf']);
end
