function [pBirth]= calcula_pBirth(X,chi,tau,nu,alpha1,gamma1,gamma2,beta,Nsim,Snew,nnew,rnew,S,n,r,md)
%% Gibbs sampler for Tapas Bar
%% Input
% X is the discrete observation matrix, size D (dimensionality)xT(number of observatios)
% Q number of states, minimum value 0
% s2x variance of the observations
% s2phi variance of the same state among the chains
% s20 variance of 
% alpha IBP concentration parameter
% gamma_q hyperparameter of am which is related to the probability of being
% in zeros-state
% beta_q dirichlet distribution parameter for jumping from one
% state(different from zero) to de rnew
% Nsim number of iterations of the Gibbs sampler
%% Output
% Snew inferred IBP matrix
% Phiest mean of the posterior of p(phi|rnew of variables)
% Mest number of inferred chains
% nnew number of jumps from one state to the rnew in each chain
%% Inicializacion
Q=3;
fk=zeros(1,3);
pBirth=0;

tvec=find(S(md,:)~=0);
[D T]=size(X);

sd= zeros(1,T);
p0 = [1/3 1/3 1/3];%p0/sum(p0);
for t=1:T
    sd(1,t) = find(mnrnd(1,p0)==1)-1;
    pBirth= pBirth+log(p0(sd(1,t)+1));
end
nd=zeros(Q,Q);
for q=1:3
    for k=1:3
        nd(q,k)=sum([0 sd(1,1:end-1)]==q-1 & sd(1,:)==k-1);
    end
end
rd=r(md);
pBirth=pBirth+log(poisspdf(r(md)-1,beta));
Nsim=1;
%% Inferencia
for it=1:Nsim
    
%% Muestreo de Z
    
    %% t=1:T
    for t=1:T 

        if t==1 %% Para t=1
            %t=1;   
            nd(sd(t)+1,sd(t+1)+1) = nd(sd(t)+1,sd(t+1)+1) -1;
            nd(1,sd(t)+1) = nd(1,sd(t)+1)-1;
            p=zeros(1,Q);

            if sd(t+1)==0
                p(1)= log(nd(1,1)+alpha1+1)+log(nd(1,1)+alpha1)-log(1+alpha1+sum(nd(1,:)))-log(alpha1+sum(nd(1,:)));
                p(2:Q)= -log(2)+log(sum(nd(1,2:end)))+log(gamma1+nd(2:end,1)')-log(gamma1+gamma2+sum(nd(2:Q,:),2).')-log(alpha1+sum(nd(1,:)));
            else
                p(1)=  -log(2)+log(alpha1+nd(1,1))-log(alpha1+sum(nd(1,:)))-log(1+alpha1+sum(nd(1,:)));
                p(2:Q)=  -log(4)+log(gamma2+sum(nd(2:end,2:end),2)')-log(gamma1+gamma2+sum(nd(2:Q,:),2).')-log(alpha1+sum(nd(1,:)));
            end
            for k=1:Q
                sd(t)=k-1;
                fk(k)= mLik_TB(X,[Snew;sd],[rnew rd],chi,tau,nu);
            end
            prob = fk+ p;
            maximo = max(prob);
            prob=exp(prob-maximo);
            prob=prob/sum(prob);
            sd(t)=S(md,t);
            pBirth=pBirth+log(prob(sd(t)+1));
            nd(sd(t)+1,sd(t+1)+1) = nd(sd(t)+1,sd(t+1)+1) +1;
            nd(1,sd(t)+1) = nd(1,sd(t)+1)+1;
                
            
        elseif t==T %% t=T

            j=sd(t-1);
            nd(j+1,sd(t)+1) = nd(j+1,sd(t)+1) -1;

            p=zeros(1,Q);
            if(j==0)
                p(1)= log(nd(1,1)+alpha1)-log(alpha1+sum(nd(1,:)));
                p(2:Q)= -log(2)+log(sum(nd(1,2:end)))-log(alpha1+sum(nd(1,:)));
            else
                p(1)= log(gamma1+nd(j+1,1))- log(gamma1+gamma2+sum(nd(j+1,:),2));
                p(2:Q) = -log(2)+log(gamma2+nd(j+1,2:Q))- log(gamma1+gamma2+sum(nd(j+1,:),2));
            end

            for k=1:Q
                sd(t)=k-1;
                fk(k)= mLik_TB(X,[Snew;sd],[rnew rd],chi,tau,nu);
            end
            prob = fk+ p;
            maximo = max(prob);
            prob=exp(prob-maximo);
            prob=prob/sum(prob);
            sd(t)=S(md,t);
            pBirth=pBirth+log(prob(sd(t)+1));
            nd(j+1,sd(t)+1) = nd(j+1,sd(t)+1) +1;
                
            
            
        else %t~=1 t~=T 
            j=sd(t-1);
            nd(j+1,sd(t)+1) = nd(j+1,sd(t)+1) -1;
            nd(sd(t)+1,sd(t+1)+1) = nd(sd(t)+1,sd(t+1)+1) -1;

            p=zeros(1,Q);
            if(j==0)
                if sd(t+1)==0
                    p(1)= log(nd(1,1)+alpha1+1)+log(nd(1,1)+alpha1)-log(alpha1+1+sum(nd(1,:)))-log(alpha1+sum(nd(1,:)));
                    p(2:Q)= -log(2)+log(sum(nd(1,2:end)))+log(gamma1+nd(2:end,1)')-log(gamma1+gamma2+sum(nd(2:Q,:),2).')-log(alpha1+sum(nd(1,:)));
                else
                    p(1)=  -log(2)+log(alpha1+nd(1,1))-log(alpha1+sum(nd(1,:)))-log(alpha1+1+sum(nd(1,:)));
                    p(2:Q)=  -log(4)+log(gamma2+sum(nd(2:end,2:end),2)')-log(gamma1+gamma2+sum(nd(2:Q,:),2).')-log(alpha1+sum(nd(1,:)));
                end

            else
                if sd(t+1)==0
                    p(1)= log(gamma1+nd(j+1,1))- log(gamma1+gamma2+sum(nd(j+1,:),2))-log(alpha1+sum(nd(1,:)))+log(alpha1+nd(1,1));
                    p(2:Q) = -log(2)+log(gamma1+nd(2:Q,1)')+log(gamma2+sum(nd(j+1,2:Q)))-log(gamma1+gamma2+sum(nd(2:Q,:),2)')-log(gamma1+gamma2+sum(nd(j+1,:))+([2:Q]==j+1));
                else
                    p(1)= -log(2)+ log(gamma1+nd(j+1,1))- log(gamma1+gamma2+sum(nd(j+1,:),2))+log (sum(nd(1,2:Q)))-log(alpha1+sum(nd(1,:)));
                    p(2:Q) = -log(4)+log(gamma2+sum(nd(2:end,2:end),2)')-log(gamma1+gamma2+sum(nd(2:Q,:),2).')+ log(gamma2+sum(nd(j+1,2:end),2)'+([2:Q]==j+1))-log(gamma1+gamma2+sum(nd(j+1,:),2).'+([2:Q]==j+1));
                end
            end
            for k=1:Q
                sd(t)=k-1;
                fk(k)= mLik_TB(X,[Snew;sd],[rnew rd],chi,tau,nu);

            end
            prob = fk+ p;
            maximo = max(prob);
            prob=exp(prob-maximo);
            prob=prob/sum(prob);
            sd(t)=S(md,t);
            pBirth=pBirth+log(prob(sd(t)+1));
            nd(j+1,sd(t)+1) = nd(j+1,sd(t)+1) +1;
            nd(sd(t)+1,sd(t+1)+1) = nd(sd(t)+1,sd(t+1)+1) +1;
               
            
        end

        
    end
end
    Zest=zeros(size(S));
    Zest([1:md-1 md+1:end],:)=Snew;
    Zest(md,:)=sd;
    Mest=size(Zest,1);
    
    nest=zeros(size(n));
    nest(:,:,[1:md-1 md+1:end])=nnew;
    nest(:,:,md)=nd;
    
    rest=zeros(size(r));
    rest([1:md-1 md+1:end])=rnew;
    rest(md)=rd;
    pBirth=log(poisspdf(r(md)-1,beta));
    %% Para el resto
     for t=tvec 

        if t==1 %% Para t=1
            for m=[1:md-1 md+1:Mest]
                %t=1;   
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) -1;
                nest(1,Zest(m,t)+1,m) = nest(1,Zest(m,t)+1,m)-1;
                p=zeros(1,Q);

                if Zest(m,t+1)==0
                    p(1)= log(nest(1,1,m)+alpha1+1)+log(nest(1,1,m)+alpha1)-log(1+alpha1+sum(nest(1,:,m)))-log(alpha1+sum(nest(1,:,m)));
                    p(2:Q)= -log(2)+log(sum(nest(1,2:end,m)))+log(gamma1+nest(2:end,1,m)')-log(gamma1+gamma2+sum(nest(2:Q,:,m),2).')-log(alpha1+sum(nest(1,:,m)));
                else
                    p(1)=  -log(2)+log(alpha1+nest(1,1,m))-log(alpha1+sum(nest(1,:,m)))-log(1+alpha1+sum(nest(1,:,m)));
                    p(2:Q)=  -log(4)+log(gamma2+sum(nest(2:end,2:end,m),2)')-log(gamma1+gamma2+sum(nest(2:Q,:,m),2).')-log(alpha1+sum(nest(1,:,m)));
                end

                for k=1:Q
                    Zest(m,t)=k-1;
                    fk(k)= mLik_TB(X,Zest,rest,chi,tau,nu);
                end
                prob = fk+ p;
                maximo = max(prob);
                prob=exp(prob-maximo);
                prob=prob/sum(prob);
                
                Zest(m,t) = S(m,t);
                pBirth=pBirth+log(prob(Zest(m,t)+1));
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) +1;
                nest(1,Zest(m,t)+1,m) = nest(1,Zest(m,t)+1,m)+1;
            end
            
        elseif t==T %% t=T
            for m=[1:md-1 md+1:Mest]
                j=Zest(m,t-1);
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) -1;

                p=zeros(1,Q);
                if(j==0)
                    p(1)= log(nest(1,1,m)+alpha1)-log(alpha1+sum(nest(1,:,m)));
                    p(2:Q)= -log(2)+log(sum(nest(1,2:end,m)))-log(alpha1+sum(nest(1,:,m)));
                else
                    p(1)= log(gamma1+nest(j+1,1,m))- log(gamma1+gamma2+sum(nest(j+1,:,m),2));
                    p(2:Q) = -log(2)+log(gamma2+nest(j+1,2:Q,m))- log(gamma1+gamma2+sum(nest(j+1,:,m),2));
                end

                for k=1:Q
                    Zest(m,t)=k-1;
                    fk(k)= mLik_TB(X,Zest,rest,chi,tau,nu);
                end
                prob = fk+ p;
                maximo = max(prob);
                prob=exp(prob-maximo);
                prob=prob/sum(prob);
                
                Zest(m,t) = S(m,t);
                pBirth=pBirth+log(prob(Zest(m,t)+1));
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) +1;
            end
                
            
            
        else %t~=1 t~=T 
            for m=[1:md-1 md+1:Mest]
                j=Zest(m,t-1);
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) -1;
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) -1;

                p=zeros(1,Q);
                if(j==0)
                    if Zest(m,t+1)==0
                        p(1)= log(nest(1,1,m)+alpha1+1)+log(nest(1,1,m)+alpha1)-log(alpha1+1+sum(nest(1,:,m)))-log(alpha1+sum(nest(1,:,m)));
                        p(2:Q)= -log(2)+log(sum(nest(1,2:end,m)))+log(gamma1+nest(2:end,1,m)')-log(gamma1+gamma2+sum(nest(2:Q,:,m),2).')-log(alpha1+sum(nest(1,:,m)));
                    else
                        p(1)=  -log(2)+log(alpha1+nest(1,1,m))-log(alpha1+sum(nest(1,:,m)))-log(alpha1+1+sum(nest(1,:,m)));
                        p(2:Q)=  -log(4)+log(gamma2+sum(nest(2:end,2:end,m),2)')-log(gamma1+gamma2+sum(nest(2:Q,:,m),2).')-log(alpha1+sum(nest(1,:,m)));
                    end

                else
                    if Zest(m,t+1)==0
                        p(1)= log(gamma1+nest(j+1,1,m))- log(gamma1+gamma2+sum(nest(j+1,:,m),2))-log(alpha1+sum(nest(1,:,m)))+log(alpha1+nest(1,1,m));
                        p(2:Q) = -log(2)+log(gamma1+nest(2:Q,1,m)')+log(gamma2+sum(nest(j+1,2:Q,m)))-log(gamma1+gamma2+sum(nest(2:Q,:,m),2)')-log(gamma1+gamma2+sum(nest(j+1,:,m))+([2:Q]==j+1));
                    else
                        p(1)= -log(2)+ log(gamma1+nest(j+1,1,m))- log(gamma1+gamma2+sum(nest(j+1,:,m),2))+log (sum(nest(1,2:Q,m)))-log(alpha1+sum(nest(1,:,m)));
                        p(2:Q) = -log(4)+log(gamma2+sum(nest(2:end,2:end,m),2)')-log(gamma1+gamma2+sum(nest(2:Q,:,m),2).')+ log(gamma2+sum(nest(j+1,2:end,m),2)'+([2:Q]==j+1))-log(gamma1+gamma2+sum(nest(j+1,:,m),2).'+([2:Q]==j+1));
                    end
                end
                for k=1:Q
                    Zest(m,t)=k-1;
                    fk(k)= mLik_TB(X,Zest,rest,chi,tau,nu);

                end
                prob = fk+ p;
                maximo = max(prob);
                prob=exp(prob-maximo);
                prob=prob/sum(prob);
                
                Zest(m,t) = S(m,t);
                pBirth=pBirth+log(prob(Zest(m,t)+1));
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) +1;
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) +1;
            end
               
            
        end

        
    end


