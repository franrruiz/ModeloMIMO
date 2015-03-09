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
for i=1:length(sweep_vec)
    hold on;
    idxGood = find(ALL_MEST(:,end,i)==Ntaux);
    idxBad = find(ALL_MEST(:,end,i)~=Ntaux);
    if(strcmp(sweep_var,'Nt'))
        idxGood = find(ALL_MEST(:,end,i)==sweep_vec(i));
        idxBad = find(ALL_MEST(:,end,i)~=sweep_vec(i));
    end
    if(~isempty(idxBad))
        plot(sweep_vec(i),ALL_MMSE(idxBad,end,i),'x','MarkerEdgeColor','r','MarkerSize',5,'MarkerFaceColor','r','LineWidth',1.5);
    end
    if(~isempty(idxGood))
        plot(sweep_vec(i),ALL_MMSE(idxGood,end,i),'x','MarkerEdgeColor',[0 0.6 0],'MarkerSize',5,'MarkerFaceColor',[0 0.6 0],'LineWidth',1.5);
    end
end


