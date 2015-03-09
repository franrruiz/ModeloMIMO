r_ = cell(2,2,1,8,100);
Mfin_ = zeros(2,2,1,8,100);
MSE_ = zeros(2,2,1,8,100);
BER_ = zeros(2,2,1,8,100);
pOnOff_ = zeros(2,2,1,8,100);
acierto = zeros(2,2,1,8,100);
aciertoM = zeros(2,2,1,8,100);

T=200;

for f=0:1
    for fChi=0:1
        for l=20
            for snr=-5:2:10
                
                BASEDIR=['/export/clusterdata/ivalera/modMIMOpruebasCluster/flag' num2str(f) '_flagChi' num2str(fChi) '_L' num2str(l) '_SNR' num2str(snr) '_T'  num2str(T)];

                if ( isdir(BASEDIR) == 0 )
                    mkdir(BASEDIR);
                end

                %if exist([BASEDIR '/../Tarea1Bis/resultadosQini' num2str(Qini) '_N' num2str(N) '_L' num2str(L) '_M' num2str(M) '_SNR' num2str(SNRdB) '_T'  num2str(T) '.mat'],'file')==0    

                for itcluster=1:100
                        load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);
                       
                        r_{f+1,fChi+1,l-19,snr*0.5+3.5,itcluster} = rest;
                        Mfin_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster) = size(Zest,1);
                        
                        [raux idxr]=sort(rest);
                        [Maux idxM]=sort(Mvec);
                        
                        if size(Zest,1)==N
                            aciertoM(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)=1;
                        end
                        if size(Zest,1)==N && sum(abs(raux-Maux))==0
                            
                            acierto(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)=1;
                            Hest=Hest';
                            Haux=zeros(L,max(rest)*size(Zest,1));
                            Hestaux=zeros(L,max(rest)*size(Zest,1));
                            for m=1:max(Mvec)
                                Haux(:,(max(Mvec)-m)*size(Zest,1)+(1:length(idxr)))=H(:,(m-1)*size(Zest,1)+idxM);
                                Hestaux(:,(max(Mvec)-m)*size(Zest,1)+(1:length(idxr)))=Hest(:,(max(Mvec)-m)*size(Zest,1)+idxr);
                            end
                            MSE_aux=zeros(1,2);
                            signo=zeros(1,size(Zest,1));
                            Zaux=zeros(size(Zest));
                            for m=1:size(Zest,1)
                                MSE_aux(1)=sum(sum((Haux(:,m:size(Zest,1):(raux(m)-1)*size(Zest,1)+m)-Hestaux(:,m:size(Zest,1):(raux(m)-1)*size(Zest,1)+m)).^2));
                                MSE_aux(2)=sum(sum((Haux(:,m:size(Zest,1):(raux(m)-1)*size(Zest,1)+m)+Hestaux(:,m:size(Zest,1):(raux(m)-1)*size(Zest,1)+m)).^2));
                                
                                [val signo(m)]=min(MSE_aux);
                                MSE_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)=MSE_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)+val;
                                if signo(m)==1
                                    Zaux(m,find(Zest(idxr(m),:)==1))=1;
                                    Zaux(m,find(Zest(idxr(m),:)==2))=-1;
                                else
                                    Zaux(m,find(Zest(idxr(m),:)==1))=-1;
                                    Zaux(m,find(Zest(idxr(m),:)==2))=1;
                                end
                            end
                                MSE_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)=MSE_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)/(size(Zest,1)*sum(rest));
                                simb=simb(idxM,max(Mvec):end);
                                idxSimb = find(simb~=0);
                                BER_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)=sum(sum(abs(simb(idxSimb)-Zest(idxSimb))))/length(idxSimb);
                                
                                pOnOff_(f+1,fChi+1,l-19,snr*0.5+3.5,itcluster)=sum(sum(abs(Zaux)~=abs(simb)))/(N*T);
                        end
                end
            end
        end
    end
end
                        
                        