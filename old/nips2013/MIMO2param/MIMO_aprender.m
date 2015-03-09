function  [Zest,Hest,MestIt,pX_ZIt,rest,nest] = MIMO_aprender(X,s2x,s2h,alpha1,alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, sim)

MestIt=zeros(1,Nsim);
pX_ZIt=zeros(1,Nsim);
Q=3;
[D T]=size(X);
%% Inicializacion

%Mest=poissrnd(alpha);
Mest=5;%(rand>0.5)+1;
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
   [Zest, rest,MestIt(it),pX_ZIt(it),nest]= MIMO_GibbsZ(X,s2x,s2h,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,1,Zest,nest,rest);
    %% Muestreo r
    for m=1:MestIt(it)
        rnew=rest;
        rnew(m)=poissrnd(rest(m)-1)+1;
        pA=mLik_TB(X,Zest,rnew,s2x,s2h)-mLik_TB(X,Zest,rest,s2x,s2h)+ log(poisspdf(rnew(m)-1,beta))-log(poisspdf(rest(m)-1,beta))+ log(poisspdf(rest(m)-1,rnew(m)-1))-log(poisspdf(rnew(m)-1,rest(m)-1));
        
        if rand<exp(pA)
            rest=rnew;
        end
    end
    
    %%
	pX_ZIt(it) = mLik_TB(X,Zest,rest,s2x,s2h);

    toc    
    MestIt(it)
    rest
    save(['Sim' num2str(sim) 'It' num2str(it)], 'rest', 'Zest')

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
W=1/s2x*ZR'*ZR+1/s2h*eye(Mest*K);

Hest=1/s2x*(W\ZR')*X';


end