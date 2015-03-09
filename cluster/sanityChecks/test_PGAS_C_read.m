AVG_PPGAS1 = zeros(2,30,5);
AVG_PPGAS2 = zeros(2,30,5);

vecItCluster = 1:10;
for itCluster=vecItCluster
    load(['/export/clusterdata/franrruiz87/ModeloMIMO/results/synthetic/pruebas/results' num2str(itCluster) '.mat'],'P_PGAS1','P_PGAS2');
    
    AVG_PPGAS1 = AVG_PPGAS1+P_PGAS1/length(vecItCluster);
    AVG_PPGAS2 = AVG_PPGAS2+P_PGAS2/length(vecItCluster);
end


disp(['Mean difference: ' num2str(mean(abs(AVG_PPGAS1(:)-AVG_PPGAS2(:))))]);
[maxDif idx] = max(abs(AVG_PPGAS1(:)-AVG_PPGAS2(:)));
disp(['Max difference: ' num2str(maxDif) ' on a probability of ' num2str(AVG_PPGAS1(idx))]);

hist((AVG_PPGAS1(:)-AVG_PPGAS2(:)),100);

[valnul idx] = sort(abs(AVG_PPGAS1(:)-AVG_PPGAS2(:)),'descend');

[a b c] = ind2sub([2 30 5],idx(1:7))
for i=1:7
    disp([AVG_PPGAS1(idx(i)) AVG_PPGAS2(idx(i))])
end

