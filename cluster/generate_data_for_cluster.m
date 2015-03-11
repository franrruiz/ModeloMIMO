clear all;

simId = 11;
folderToSave = ['/export/clusterdata/franrruiz87/ModeloMIMO/data/syn/' num2str(simId) '/'];

Nt = 50;
Nr = 200;
L = 20;
Mmax = 10;

for T=[500 1000]
    for i=1:100
        fileToSave = [folderToSave 'T' num2str(T) '_itCluster' num2str(i) '.mat'];
        symbolsGen = cell(1,Mmax);
        seqGen = cell(1,Mmax);

        for M=1:Mmax
            constellation = qammod(0:2^M-1,2^M,[],'gray');
            constellation = constellation/sqrt(mean(abs(constellation.^2)));
            auxConstellation = [0 constellation];
            seqGen{M} = randint(Nt,T,[1 length(constellation)]);
            for nt=1:Nt
                Tini = randint(1,1,[1 round(T/2)]);
                Tend = min(Tini+round(T/2)-1,T);
                seqGen{M}(nt,setdiff(1:T,Tini:Tend)) = 0;
            end
            symbolsGen{M} = zeros(size(seqGen{M}));
            symbolsGen{M}(:) = auxConstellation(1+seqGen{M});
        end
        
        channelGen = randn(Nr,Nt,L)+1i*randn(Nr,Nt,L);
        
        noiseGen = randn(Nr,T)+1i*randn(Nr,T);
        
        save(fileToSave,'symbolsGen','seqGen','channelGen','noiseGen');
    end
end



