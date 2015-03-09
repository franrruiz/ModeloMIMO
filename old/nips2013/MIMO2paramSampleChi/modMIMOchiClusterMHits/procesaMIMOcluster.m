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
mem_ = zeros(1,1,100,N);

flag=0;
BASEDIR=['/export/clusterdata/ivalera/modMIMOsampleChiNewMHitsN' num2str(N) 'T' num2str(T) '/flag' num2str(flag) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T)];
if ( isdir(BASEDIR) == 0 )
    mkdir(BASEDIR);
end

if exist([BASEDIR '/../ResultadosBER/N' num2str(N) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T) '.mat'],'file')==0
    for itcluster=1:100
        load([BASEDIR '/SimFinal_' num2str(itcluster) '.mat']);

        idxCut=find(sum(Zest~=0,2)>=40);
        Zest = Zest(idxCut,:);
        rest=rest(idxCut);
        Zest(find(Zest==2))=-1;
        
        
        Mest=length(idxCut);  

        Nfin_(itcluster) = length(rest);
        if length(rest)==N
            Mfin_(itcluster) = sum(sort(rest)==sort(Mvec))==N;
        end
        
        simb=simb(:,max(Mvec):end);
        Mt=N;
        N=max(Mt,Mest);
        if Mt>Mest
            Zest=[Zest; zeros(Mt-Mest,T)];
            rest=[rest zeros(1,Mt-Mest)];
        elseif Mt<Mest
            simb=[simb; zeros(Mest-Mt,T)];
            Mvec=[Mvec zeros(1,Mest-Mt)];
        end
        
        tabOrdenTx = perms(1:N);
        combBERs = Inf;
        Zaux = Zest;
        raux=rest;
       
        for ordenTx=1:factorial(N)
            Zord =Zest(tabOrdenTx(ordenTx,:),:);   	
            for i=0:2^N-1
                signos=de2bi(i,N);
                signos(find(signos==0))=-1;
                
                errBERs=sum(sum(simb~=Zord.*repmat(signos',1,T)));
                if errBERs< combBERs
                     combBERs=errBERs;
                     Zaux=Zord.*repmat(signos',1,T);
                     raux=rest(tabOrdenTx(ordenTx,:));
                end
            end
        end 
        
        %% Calcula la media de las H's
        Zest=Zaux;
        Mest=N;
        rest=raux;
        K=max(max(rest),max(Mvec));
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
        
        rmax=zeros(1,N);
        HestTx=cell(1,N);
        HTx=cell(1,N);
        errMMSE=0;
        for n=1:Mt
            %HestTx{n} = zeros(L,K);
            HTx{n}=fliplr(H(:,n:Mt:end));
            HTx{n}=[HTx{n} zeros(D,K-max(Mvec))];
        end
        for n=Mt+1:N
            HTx{n}=zeros(D,K);
        end
        for n=1:N
            %HestTx{n} = zeros(L,K);
            %HTx{n}=fliplr(H(:,n:Mt:end));
            %HTx{n}=[HTx{n} zeros(D,K-max(Mvec))];
            HestTx{n}=Hest(n:N:end,:)';
            errMMSE=errMMSE+sum(sum((HTx{n}(:,1:max(Mvec(n),raux(n)))-HestTx{n}(:,1:max(Mvec(n),raux(n)))).^2));    
            rmax(n)=max(Mvec(n),raux(n)); 
        end
        
        if(Mt==Nfin_(itcluster))
            mem_(1,1,itcluster,:) = Mvec-raux;
        end
        
        BER_(itcluster) = sum(sum(Zaux~=simb))/(N*T);%sum(Zaux(idxSimb)~=simb(idxSimb))/length(idxSimb);
        MSE_(itcluster) = errMMSE/sum(rmax);
        pOnOff_(itcluster) = sum(sum(abs(Zaux)~=abs(simb)))/(N*T);

    end

    save ([BASEDIR '/../ResultadosBER/N' num2str(Mt) '_L' num2str(L) '_SNR' num2str(SNR) '_T'  num2str(T) '.mat'], 'BER_','MSE_','Mfin_','Nfin_','pOnOff_','mem_')
end

    