function procesaMIMOcluster(N,L,SNR,T)

addpath '/export/clusterdata/ivalera/minFunc_2012/'
addpath '/export/clusterdata/ivalera/minFunc_2012/autoDif'
addpath '/export/clusterdata/ivalera/minFunc_2012/logisticExample'
addpath '/export/clusterdata/ivalera/minFunc_2012/minFunc'
addpath '/export/clusterdata/ivalera/minFunc_2012/minFunc/mex'
addpath '/export/clusterdata/ivalera/minFunc_2012/minFunc/compiled'

%optionsMin = optimset('GradObj','on','Display','off','TolX',1e-5,'TolFun',1e-10);
options = [];
options.display = 'none';
options.useMex=0;
options.Method = 'lbfgs';
options.MaxIter=1000;
options.progTol=1e-10;



Nfin_ = zeros(1,1,100);%numero de transmisores estimado
Mfin_ = zeros(1,1,100);%binario
MSE_ = zeros(1,1,100);
BER_ = zeros(1,1,100);
pOnOff_ = zeros(1,1,100);
flag=0;
BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(flag) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
if ( isdir(BASEDIR) == 0 )
    mkdir(BASEDIR);
end

if exist([BASEDIR '/../ResultadosNew/N' num2str(N) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T) '.mat'],'file')==0
    for itcluster=1:100
        load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);

        idxCut=find(sum(Zest~=0,2)>=40);
        Zest = Zest(idxCut,:);
        rest=rest(idxCut);
        Zest(find(Zest==2))=-1;
        
        %% Calcula la media de las H's
        Mest=length(idxCut);
        K=max(rest);
        S=Zest';
        %S(find(S==2))=-1;
        Z=zeros(T,Mest*K);
        R=zeros(Mest*K,Mest*K);
        for k=1:K
            Z(k:end,(k-1)*Mest+1:k*Mest)=S(1:T-k+1,:);
            R((k-1)*Mest+1:k*Mest,(k-1)*Mest+1:k*Mest)=diag(rest>=k);
        end
        ZR=Z*R;
        Hest = zeros(Mest*K,D);
        %Hest = fminunc(@(H) posteriorPhiPdf_H_s2x(H,X',ZR,chi,tau,nu), 10*randn(size(Hest)), optionsMin);
        Hest = minFunc(@(H) posteriorPhiPdf_H_s2x_2(H,X',ZR,chi,tau,nu), 10*randn(size(Hest(:))),options);
        Hest=reshape(Hest,Mest*K,D);

        Nfin_(itcluster) = length(rest);
        if length(rest)==N
            Mfin_(itcluster) = sum(sort(rest)==sort(Mvec))==N;
        end

        if length(rest)==N && max(rest)==max(Mvec)         
            tabOrdenTx = perms(1:N);
            combMMSEs = Inf;
            Zaux = Zest;
            raux=rest;
            HestTx=cell(1,N);
            Haux=cell(1,N);
            HTx=cell(1,N);
            for n=1:N
                %HestTx{n} = zeros(L,K);
                HTx{n}=fliplr(H(:,n:N:end));
            end
            for ordenTx=1:factorial(N)
                for n=1:N
                    Haux{n} =Hest(tabOrdenTx(ordenTx,n):N:end,:)';   	
                end
                for i=0:2^N-1
                    errH=0;
                    signos=de2bi(i,N);
                    signos(find(signos==0))=-1;
                    for n=1:N
                        errH=errH+sum(sum(abs(HTx{n}-signos(n)*Haux{n}).^2));
                    end
                    if errH< combMMSEs
                        combMMSEs=errH;
                        for n=1:N
                            HestTx{n}=Haux{n}*signos(n);  
                            Zaux(n,:)=Zest(tabOrdenTx(ordenTx,n),:)*signos(n);
                            raux(n)=rest(tabOrdenTx(ordenTx,n));
                        end
                        
                    end
                end
            end 


            %combMMSEs %sum(Zaux~=simb(:,3:end),2)
            rmax=zeros(1,N);
            for n=1:N
               rmax(n)=max(Mvec(n),raux(n)); 
            end
            simb=simb(:,K:end);
            idxSimb = find(simb~=0);
            BER_(itcluster) = sum(sum(Zaux~=simb))/(N*T);%sum(Zaux(idxSimb)~=simb(idxSimb))/length(idxSimb);
            MSE_(itcluster) = combMMSEs/sum(rmax);
            pOnOff_(itcluster) = sum(sum(abs(Zaux)~=abs(simb)))/(N*T);
        end
    end

    save ([BASEDIR '/../ResultadosNew/N' num2str(N) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T) '.mat'], 'BER_','MSE_','Mfin_','Nfin_','pOnOff_')
end

    