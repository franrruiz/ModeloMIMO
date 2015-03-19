for SNR=-18:3:0
    SNR
    for itCluster = 1:50
        addIndividualResultsToSim(1000,5,20,2,1,1,SNR,10000,0,0,itCluster,11);
    end
end

for Nr=5:5:25
    Nr
    for itCluster = 1:50
        addIndividualResultsToSim(1000,5,Nr,2,1,1,-3,10000,0,0,itCluster,11);
    end
end

for M=2:7
    M
    for itCluster = 1:50
        addIndividualResultsToSim(1000,5,20,M,1,1,-3,10000,0,0,itCluster,11);
    end
end

for Nt=2:2:10
    Nt
    for itCluster = 1:50
        addIndividualResultsToSim(1000,Nt,20,2,1,1,-3,10000,0,0,itCluster,11);
    end
end






