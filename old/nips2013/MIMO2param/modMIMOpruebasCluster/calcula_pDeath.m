function [pDeath]= calcula_pDeath(X,chi,tau,nu,alpha1,gamma1,gamma2,Nsim,Snew,nnew,rnew,S,n,r)
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
% tvec
Nsim=1;
[D T]=size(X);
Mest=size(Snew,1);
%% Inicializacion
Q=3;
fk=zeros(1,3);
pDeath=-log(Mest);
tvec=find(Snew(end,:)~=0);

Zest=Snew(1:end-1,:);
nest=nnew(:,:,1:end-1);
rest=rnew(1:end-1);
%% Inferencia
for it=1:Nsim
    
%% Muestreo de Z
    
    %% t=1:T
    for t=tvec  

        if t==1 %% Para t=1
           
            for m=1:Mest-1
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
                Zest(m,t) = S(m,t);
                pDeath=pDeath+log(prob(Zest(m,t)+1));
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) +1;
                nest(1,Zest(m,t)+1,m) = nest(1,Zest(m,t)+1,m)+1;
                
            end
        elseif t==T %% t=T
           
            for m=1:Mest-1
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
                pDeath=pDeath+log(prob(Zest(m,t)+1));
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) +1;
                
            end
            
        else %t~=1 t~=T 
            
            for m=1:Mest-1
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
                pDeath=pDeath+log(prob(Zest(m,t)+1));
                nest(j+1,Zest(m,t)+1,m) = nest(j+1,Zest(m,t)+1,m) +1;
                nest(Zest(m,t)+1,Zest(m,t+1)+1,m) = nest(Zest(m,t)+1,Zest(m,t+1)+1,m) +1;
               
            end
        end

        
    end
end



