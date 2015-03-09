function  [Zest,Hest,MestIt,pX_ZIt,chiIt,rest,nest] = MIMO_aprender_H_s2x(X,chi,tau,nu,alpha1,alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, sim)


MestIt=zeros(1,Nsim);
pX_ZIt=zeros(1,Nsim);
chiIt=zeros(1,Nsim);
Q=3;
[D T]=size(X);

addpath '/Users/ivalera/Desktop/Tesis Doctoral/codigo/minFunc_2012/'
addpath '/Users/ivalera/Desktop/Tesis Doctoral/codigo/minFunc_2012/autoDif'
addpath '/Users/ivalera/Desktop/Tesis Doctoral/codigo/minFunc_2012/logisticExample'
addpath '/Users/ivalera/Desktop/Tesis Doctoral/codigo/minFunc_2012/minFunc'
addpath '/Users/ivalera/Desktop/Tesis Doctoral/codigo/minFunc_2012/minFunc/mex'
addpath '/Users/ivalera/Desktop/Tesis Doctoral/codigo/minFunc_2012/minFunc/compiled'

%optionsMin = optimset('GradObj','on','Display','off','TolX',1e-5,'TolFun',1e-10);
options = [];
options.display = 'none';
options.useMex=0;
options.Method = 'lbfgs';
options.MaxIter=1000;
options.progTol=1e-10;
%% Inicializacion

%Mest=poissrnd(alpha);
Mest=1;%(rand>0.5)+1;
Zest=zeros(Mest,T);
nest = zeros(Q,Q,Mest);


%p0=ones(1,Q)/(Q);
%p0 = ones(1,Q)+[1 0.5 zeros(1,Q-2)];
p0 = [0.8 0.1 0.1];%p0/sum(p0);

for m=1:Mest
    Zest(m,1)=find(mnrnd(1,p0)==1)-1;
    nest(1,Zest(m,1)+1,m) = 1;
    for t=2:T
        Zest(m,t) = find(mnrnd(1,p0)==1)-1;
        nest(Zest(m,t-1)+1,Zest(m,t)+1,m) =  nest(Zest(m,t-1)+1,Zest(m,t)+1,m) + 1;
    end
end
rest=poissrnd(beta,1,Mest)+1;
% nest=n;
% Mest=M;
% Zest=Z;
% thetaest=theta;

for it=1:Nsim
    it
    tic
    %% pasos a) y b): Gibbs
   [Zest, rest,MestIt(it),pX_ZIt(it),nest]= MIMO_GibbsZ_H_s2x(X,chi,tau,nu,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,1,Zest,nest,rest);
    %% Muestreo r
    for m=1:MestIt(it)
        rnew=rest;
        rnew(m)=poissrnd(rest(m)/2)+1;
        if rand<0.5 %izqda (sin deplazamiento)
            pA=mLik_TB(X,Zest,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu)+ log(poisspdf(rnew(m)-1,beta))-log(poisspdf(rest(m)-1,beta))+ log(poisspdf(rest(m)-1,rnew(m)/2))-log(poisspdf(rnew(m)-1,rest(m)/2));
            Znew=Zest;
            nnew=nest;
        else
            Znew=Zest;
            if rnew(m)>rest(m)
                Znew(m,end-abs(rnew(m)-rest(m))+1:end)=(mnrnd(1,[1/3 1/3 1/3],abs(rnew(m)-rest(m)))*[0 1 2]')';
                Znew(m,1:end-abs(rnew(m)-rest(m)))= Zest(m,abs(rnew(m)-rest(m))+1:end);
                
            else
                Znew(m,1:abs(rnew(m)-rest(m)))=(mnrnd(1,[1/3 1/3 1/3],abs(rnew(m)-rest(m)))*[0 1 2]')';
                Znew(m,abs(rnew(m)-rest(m))+1:end)= Zest(m,1:end-abs(rnew(m)-rest(m)));
            end
            
            nnew=zeros(Q,Q,MestIt(it));
            for mm=1:MestIt(it)
                for q=1:Q
                    for k=1:Q
                        nnew(q,k,mm)=sum([0 Znew(mm,1:end-1)]==q-1 & Znew(mm,:)==k-1);
                    end
                end
            end

            pA=mLik_TB(X,Znew,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu)+ log(poisspdf(rnew(m)-1,beta))-log(poisspdf(rest(m)-1,beta))+ calcula_pS(nnew,alpha1,alpha2,gamma1,gamma2,T)-calcula_pS(nest,alpha1,alpha2,gamma1,gamma2,T)+log(poisspdf(rest(m)-1,rnew(m)/2))-log(poisspdf(rnew(m)-1,rest(m)/2));
        end
     
        
        if rand<exp(pA)
            rest=rnew;
            nest=nnew;
            Zest=Znew;
        end
       
    end
    
    %% Split/Merge cadenas
    R=-Inf;
    flagSplit=-1;
   
    m1=mnrnd(1,1/MestIt(it)*ones(1,MestIt(it)))*(1:MestIt(it))';
    m2=mnrnd(1,1/MestIt(it)*ones(1,MestIt(it)))*(1:MestIt(it))';
    if m1==m2
        flagSplit=1;
        m=m1;
    else
        flagSplit=0;
    end
    pm=-2*log(MestIt(it));
    if flagSplit==0 
        [Snew nnew pMerge flagSigno]=mergeChains(Zest,nest,m1,m2);
        rnew=[rest([1:min(m1,m2)-1 min(m1,m2)+1:max(m1,m2)-1 max(m1,m2)+1:end]) max(rest(m1),rest(m2))];
        pSplit=calcula_pSplit(Zest,Snew,X,m1,m2,rest,chi,tau,nu,alpha1,gamma1,gamma2, flagSigno);
        pmnew=-2*log(MestIt(it)+1);        
        R= mLik_TB(X,Snew,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu) - calcula_pS(nest,alpha1,alpha2,gamma1,gamma2,T) + calcula_pS(nnew,alpha1,alpha2,gamma1,gamma2,T)-log(poisspdf(min(rest(m1),rest(m2))-1,beta)) +...
            pSplit-pMerge + pmnew - pm;

    elseif flagSplit==1
        [Snew nnew rnew pSplit]=splitChains(Zest,X,nest,m,rest,chi,tau,nu, alpha1,gamma1,gamma2);
        pMerge=calcula_pMerge(Snew,m);
        pmnew=-2*log(MestIt(it)-1);    

        R= mLik_TB(X,Snew,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu) - calcula_pS(nest,alpha1,alpha2,gamma1,gamma2,T) + calcula_pS(nnew,alpha1,alpha2,gamma1,gamma2,T)+log(poisspdf(min(rnew(m),rnew(end))-1,beta)) +...
           pMerge- pSplit + pmnew - pm;
    end
    
    if rand<exp(R)
        rest=rnew;
        nest=nnew;
        Zest=Snew;
        MestIt(it)=size(Zest,1);
        if flagSplit==1
            disp('Split');
        else
            disp('Merge');
        end
    end

%     %% Birth/Death cadenas
%     R=-Inf;
%     if  MestIt(it)>1
%         flagBirth=round(rand);
%     else
%         flagBirth=1;
%     end
% 
% 
%     if flagBirth==0 
%         md=mnrnd(1,1/MestIt(it)*ones(1,MestIt(it)))*(1:MestIt(it))';
%         tvec=find(Zest(md,:)~=0);
%         [Snew, rnew, nnew, pDeath]= Gibbs_Death(X,chi,tau,nu,alpha1,gamma1,gamma2,Nsim,Zest([1:md-1 md+1:end],:),nest(:,:,[1:md-1 md+1:end]),rest([1:md-1 md+1:end]), tvec);  
%         pBirth= calcula_pBirth(X,chi,tau,nu,alpha1,gamma1,gamma2,beta,Nsim,Snew,nnew,rnew,Zest,nest,rest,md);        
%         R= mLik_TB(X,Snew,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu) - calcula_pS(nest,alpha1,alpha2,gamma1,gamma2,T) + calcula_pS(nnew,alpha1,alpha2,gamma1,gamma2,T)-log(poisspdf(rest(md)-1,beta))+...
%            pBirth-pDeath;
% 
%     elseif flagBirth==1
%         [Snew, rnew, nnew, pBirth]= Gibbs_Birth(X,chi,tau,nu,alpha1,gamma1,gamma2,beta,Nsim,Zest,nest,rest);
%         pDeath= calcula_pDeath(X,chi,tau,nu,alpha1,gamma1,gamma2,Nsim,Snew,nnew,rnew,Zest,nest,rest); 
% 
%         R= mLik_TB(X,Snew,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu) - calcula_pS(nest,alpha1,alpha2,gamma1,gamma2,T) + calcula_pS(nnew,alpha1,alpha2,gamma1,gamma2,T)+log(poisspdf(rnew(end)-1,beta)) +...
%            pDeath- pBirth;
%     end
% 
%     if rand<exp(R)
%         rest=rnew;
%         nest=nnew;
%         Zest=Snew;
%         MestIt(it)=size(Zest,1);
%         if flagBirth==1
%             disp('Birth');
%         else
%             disp('Death');
%         end
%     end
%     
    
    

    %% Sample chi
    pchi1=10;
    %pchi2=chi/pchi1;
    chinew = gamrnd(pchi1,chi/pchi1);
    R= mLik_TB(X,Zest,rest,chinew,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu)+log(gampdf(chinew,3,2))-log(gampdf(chi,3,2))+log(gampdf(chi,pchi1,chinew/pchi1))-log(gampdf(chinew,pchi1,chi/pchi1));
    %R= mLik_TB(X,Zest,rest,chinew,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu);%+log(gampdf(chinew,2,2))-log(gampdf(chi,2,2))+log(gampdf(chi,3,chinew/3))-log(gampdf(chinew,3,chi/3));
    if rand<exp(R)
        chi=chinew
    end
    %% 
	pX_ZIt(it) = mLik_TB(X,Zest,rest,chi,tau,nu);
    chiIt(it) = chi;

    toc    
    MestIt(it)
    rest
    %save(['Sim' num2str(sim) 'It' num2str(it)], 'rest', 'Zest')

end

%% Calcula la media de las H's
Mest=MestIt(it);
K=max(rest);
S=Zest';
S(find(S==2))=-1;
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

end