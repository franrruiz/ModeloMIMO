function [Snew nnew rnew pSplit]=splitChains(S,X,n,m,r,chi,tau,nu, alpha1,gamma1,gamma2)
[M,T]=size(S);
Snew=[S; zeros(1,T)];
Snew(m,:)=0;
nnew=cat(3,n,zeros(3,3));
nnew(:,:,m)=0;
m1=m;
m2=M+1;
vectorp=2*(1:r(m))/(r(m)*(r(m)+1));
rnew=[r mnrnd(1, vectorp)*(1:r(m))'];
pSplit=log(vectorp(rnew(end)));

idxN0=find(S(m,:)~=0);


t=1;
if S(m,t)==0
    nnew(1,1,m1)=nnew(1,1,m1)+1;
    nnew(1,1,m2)=nnew(1,1,m2)+1;
else
    %mLik_TB(X,S,r,chi,tau,nu)        
    p=zeros(1,5);

    s1=[0 S(m,t)*ones(1,3) 3-S(m,t)];
    s2=[S(m,t) 0 S(m,t) 3-S(m,t) S(m,t)];
    for i=1:5
        Saux=Snew;
        Saux(m1,t)=s1(i);
        Saux(m2,t)=s2(i);
        p(i)=p(i)+mLik_TB(X(:,1:t),Saux(:,1:t),rnew,chi,tau,nu);
    end
    maximo = max(p);
    p=exp(p-maximo);
    p=p/sum(p);
    j = find(mnrnd(1,p)==1);
    Snew(m1,t)=s1(j);
    Snew(m2,t)=s2(j);
    nnew(1,Snew(m1,t)+1,m1)=nnew(1,Snew(m1,t)+1,m1)+1;
    nnew(1,Snew(m2,t)+1,m2)=nnew(1,Snew(m2,t)+1,m2)+1;

    pSplit=pSplit+log(p(j));
end

for t=2:T
   
    if S(m,t)==0
        nnew(Snew(m1,t-1)+1,1,m1)=nnew(Snew(m1,t-1)+1,1,m1)+1;
        nnew(Snew(m2,t-1)+1,1,m2)=nnew(Snew(m2,t-1)+1,1,m2)+1;
    else
        %mLik_TB(X,S,r,chi,tau,nu)
        if sum(Snew(m1,1:t-1))==0            
            p=zeros(1,5);

        else
            if Snew(m1,t-1)==0
                p=[log(nnew(1,1,m1)+alpha1)-log(sum(nnew(1,:,m1))+alpha1) (log(sum(nnew(1,2:end,m1)))-log(sum(nnew(1,:,m1))+alpha1)-log(2))*ones(1,4)];

            else
                p=[log(nnew(Snew(m1,t-1)+1,1,m1)+gamma1)-log(sum(nnew(Snew(m1,t-1)+1,:,m1))+gamma1+gamma2) (log(sum(nnew(Snew(m1,t-1)+1,2:3,m1))+gamma2)-log(sum(nnew(Snew(m1,t-1)+1,:,m1))+gamma1+gamma2)-log(2))*ones(1,4)]; 

            end
        end

        if sum(Snew(m2,1:t-1))==0            
            p=p+zeros(1,5);

        else
            if Snew(m2,t-1)==0
                p=p+[(log(sum(nnew(1,2:end,m2)))-log(sum(nnew(1,:,m2))+alpha1)-log(2)) log(nnew(1,1,m2)+alpha1)-log(sum(nnew(1,:,m2))+alpha1) (log(sum(nnew(1,2:end,m2)))-log(sum(nnew(1,:,m2))+alpha1)-log(2))*ones(1,3)];

            else
                p=p+[(log(sum(nnew(Snew(m2,t-1)+1,2:3,m2))+gamma2)-log(sum(nnew(Snew(m2,t-1)+1,:,m2))+gamma1+gamma2)-log(2)) log(nnew(Snew(m2,t-1)+1,1,m2)+gamma1)-log(sum(nnew(Snew(m2,t-1)+1,:,m2))+gamma1+gamma2) (log(sum(nnew(Snew(m2,t-1)+1,2:3,m2))+gamma2)-log(sum(nnew(Snew(m2,t-1)+1,:,m2))+gamma1+gamma2)-log(2))*ones(1,3)]; 

            end
        end
        s1=[0 S(m,t)*ones(1,3) 3-S(m,t)];
        s2=[S(m,t) 0 S(m,t) 3-S(m,t) S(m,t)];
        for i=1:5
            Saux=Snew;
            Saux(m1,t)=s1(i);
            Saux(m2,t)=s2(i);
            p(i)=p(i)+mLik_TB(X(:,1:t),Saux(:,1:t),rnew,chi,tau,nu);
        end
        maximo = max(p);
        p=exp(p-maximo);
        p=p/sum(p);
        j = find(mnrnd(1,p)==1);
        Snew(m1,t)=s1(j);
        Snew(m2,t)=s2(j);
        nnew(Snew(m1,t-1)+1,Snew(m1,t)+1,m1)=nnew(Snew(m1,t-1)+1,Snew(m1,t)+1,m1)+1;
        nnew(Snew(m2,t-1)+1,Snew(m2,t)+1,m2)=nnew(Snew(m2,t-1)+1,Snew(m2,t)+1,m2)+1;

        pSplit=pSplit+log(p(j));
       
    end
end