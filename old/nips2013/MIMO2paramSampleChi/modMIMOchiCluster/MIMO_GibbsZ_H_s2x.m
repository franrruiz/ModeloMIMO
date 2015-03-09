function [Zest, rest,Mest,pX_Z,nest]= MIMO_GibbsZ_H_s2x(X,chi,tau,nu,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,Nsim,Zest,nest,rest)
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
% state(different from zero) to de rest
% Nsim number of iterations of the Gibbs sampler
%% Output
% Zest inferred IBP matrix
% Phiest mean of the posterior of p(phi|rest of variables)
% Mest number of inferred chains
% nest number of jumps from one state to the rest in each chain

[D T]=size(X);
Mest=size(Zest,1);
%% Inicializacion
Q=3;
fk=zeros(1,3);
%% Inferencia
for it=1:Nsim
    
%% Muestreo de Z
    
    %% t=1:T
    for t=1:T  

        if t==1 %% Para t=1
            m=1;
            while m<=Mest
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
                %p(2:Q) = log(beta_q+nest(2:Q,Zest(m,t+1)+1,m)')- log(Q*beta_q+sum(nest(2:Q,:,m),2)');%                 maximo = max(p);
%                 p=exp(p-maximo);
%                 p=p/sum(p);
                

                for k=1:Q
                    Zest(m,t)=k-1;
                    fk(k)= mLik_TB(X,Zest,rest,chi,tau,nu);
                end
                prob = fk+ p;
                maximo = max(prob);
                prob=exp(prob-maximo);
                prob=prob/sum(prob);
                Zest(m,t) = find(mnrnd(1,prob)==1)-1;
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) +1;
                nest(1,Zest(m,t)+1,m) = nest(1,Zest(m,t)+1,m)+1;
                %% Eliminamos cadenas vacias
                if sum(Zest(m,[1:t-1 t+1:end]))==0 && Mest>1
                    Zest=Zest([1:m-1 m+1:end],:);
                    nest= nest(:,:,[1:m-1 m+1:end]);
                    Mest=Mest-1;
                    rest=rest([1:m-1 m+1:end]);
                else             
                    m=m+1;
                end
            end
        elseif t==T %% t=T
            m=1;
            while m<=Mest
                j=Zest(m,t-1);
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) -1;
                
                p=zeros(1,Q);
                if(j==0)
                    p(1)= log(nest(1,1,m)+alpha1)-log(alpha1+sum(nest(1,:,m)));
                    p(2:Q)= -log(2)+log(sum(nest(1,2:end,m)))-log(alpha1+sum(nest(1,:,m)));
%                     maximo = max(p);
%                     p=exp(p-maximo);
%                     p=p/sum(p);
                else
                    p(1)= log(gamma1+nest(j+1,1,m))- log(gamma1+gamma2+sum(nest(j+1,:,m),2));
                    p(2:Q) = -log(2)+log(gamma2+nest(j+1,2:Q,m))- log(gamma1+gamma2+sum(nest(j+1,:,m),2));
%                     maximo = max(p);
%                     p=exp(p-maximo);
%                     p=p/sum(p);
                end

                for k=1:Q
                    Zest(m,t)=k-1;
                    fk(k)= mLik_TB(X,Zest,rest,chi,tau,nu);
                end
                prob = fk+ p;
                maximo = max(prob);
                prob=exp(prob-maximo);
                prob=prob/sum(prob);
                Zest(m,t) = find(mnrnd(1,prob)==1)-1;

                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) +1;
                %% Eliminamos cadenas vacias
                if sum(Zest(m,[1:t-1 t+1:end]))==0 && Mest>1
                    Zest=Zest([1:m-1 m+1:end],:);
                    nest= nest(:,:,[1:m-1 m+1:end]);
                    Mest=Mest-1;
                    rest=rest([1:m-1 m+1:end]);
                else           
                    m=m+1;                   
                end
            end
            
        else %t~=1 t~=T 
            m=1;
            while m<=Mest
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
                Zest(m,t) = find(mnrnd(1,prob)==1)-1;

                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) +1;
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) +1;
                %% Eliminamos cadenas vacias
                if sum(Zest(m,[1:t-1 t+1:end]))==0 && Mest>1
                    Zest=Zest([1:m-1 m+1:end],:);
                    nest= nest(:,:,[1:m-1 m+1:end]);
                    Mest=Mest-1;
                    rest=rest([1:m-1 m+1:end]);
                else         
                    m=m+1;
                end
            end
        end
        %% A?adimos cadenas
        aux=(rand<pii);
        Mnew=(1-aux)*poissrnd(alpha2*lambda/(T-1))+aux;
        rnew=[rest poissrnd(beta,1,Mnew)+1];
        snew=round(rand(Mnew,1))+1;
        Znew= [Zest; zeros(Mnew,T)];
        Znew(Mest+1:Mest+Mnew,t)=snew;
        
        pA= mLik_TB(X,Znew,rnew,chi,tau,nu)-mLik_TB(X,Zest,rest,chi,tau,nu)+Mnew*log(alpha1*alpha2)-log(factorial(sum(snew==1)))-log(factorial(sum(snew==2)))+Mnew*log(gamma1)-Mnew*log(T-1)-Mnew*log(gamma1+gamma2)-log((1-pii)*poisspdf(Mnew,alpha2*lambda/(T-1))+pii*(Mnew==1))+log((1-pii)*poisspdf(0,alpha2*lambda/(T-1)));
        acepto=rand<exp(pA);
        if acepto
            Zest=Znew;
            rest=rnew;
            nest = cat(3,nest,zeros(Q,Q,Mnew));
            
            for it_m=Mest+1:Mest+Mnew
                nest(1,Zest(it_m,t)+1,it_m) = 1;
                if t==T
                    nest(1,1,it_m) = T-1;
                else
                    nest(Zest(it_m,t)+1,1,it_m) = 1;
                    nest(1,1,it_m) = T-2;
                end
            end
            Mest=Mest+Mnew;
        end
        
    end
end

pX_Z= mLik_TB(X,Zest,rest,chi,tau,nu);
%Mest

