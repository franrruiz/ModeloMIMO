function pSplit=calcula_pSplit(S,Snew,X,m1,m2,r,chi,tau,nu,alpha1,gamma1,gamma2,flagSigno)
[M,T]=size(Snew);

n=zeros(3,3,2);
pSplit=log(2*min(r(m1),r(m2))/(max(r(m1),r(m2))*(max(r(m1),r(m2))+1)));

if flagSigno
    S(m2,find(S(m2,:)~=0))= 3-S(m2,find(S(m2,:)~=0));
end

t=1;

if Snew(end,t)==0
    n(1,1,1)=n(1,1,1)+1;
    n(1,1,2)=n(1,1,2)+1;    
else      
    p=zeros(1,5);

    s1=[0 Snew(end,t)*ones(1,3) 3-Snew(end,t)];
    s2=[Snew(end,t) 0 Snew(end,t) 3-Snew(end,t) Snew(end,t)];
    for i=1:5
        Saux=S;
        Saux(m1,t)=s1(i);
        Saux(m2,t)=s2(i);
        p(i)=p(i)+mLik_TB(X(:,1:t),Saux(:,1:t),r,chi,tau,nu);
    end
    maximo = max(p);
    p=exp(p-maximo);
    p=p/sum(p);

    pSplit=pSplit+ log(p(find(s1==S(m1,t) & s2==S(m2,t))));

    n(1,S(m1,t)+1,1)=n(1,S(m1,t)+1,1)+1;
    n(1,S(m2,t)+1,2)=n(1,S(m2,t)+1,2)+1;

end


for t=2:T  
    if Snew(end,t)==0
        n(S(m1,t-1)+1,1,1)=n(S(m1,t-1)+1,1,1)+1;
        n(S(m2,t-1)+1,1,2)=n(S(m2,t-1)+1,1,2)+1;
    else
        %mLik_TB(X,S,r,chi,tau,nu)
        if sum(S(m1,1:t-1))==0            
            p=zeros(1,5);

        else
            if S(m1,t-1)==0
                p=[log(n(1,1,1)+alpha1)-log(sum(n(1,:,1))+alpha1) (log(sum(n(1,2:end,1)))-log(sum(n(1,:,1))+alpha1)-log(2))*ones(1,4)];

            else
                p=[log(n(S(m1,t-1)+1,1,1)+gamma1)-log(sum(n(S(m1,t-1)+1,:,1))+gamma1+gamma2) (log(sum(n(S(m1,t-1)+1,2:3,1))+gamma2)-log(sum(n(S(m1,t-1)+1,:,1))+gamma1+gamma2)-log(2))*ones(1,4)]; 

            end
        end

        if sum(S(m2,1:t-1))==0            
            p=p+zeros(1,5);

        else
            if S(m2,t-1)==0
                p=p+[(log(sum(n(1,2:end,2)))-log(sum(n(1,:,2))+alpha1)-log(2)) log(n(1,1,2)+alpha1)-log(sum(n(1,:,2))+alpha1) (log(sum(n(1,2:end,2)))-log(sum(n(1,:,2))+alpha1)-log(2))*ones(1,3)];

            else
                p=p+[(log(sum(n(S(m2,t-1)+1,2:3,2))+gamma2)-log(sum(n(S(m2,t-1)+1,:,2))+gamma1+gamma2)-log(2)) log(n(S(m2,t-1)+1,1,2)+gamma1)-log(sum(n(S(m2,t-1)+1,:,2))+gamma1+gamma2) (log(sum(n(S(m2,t-1)+1,2:3,2))+gamma2)-log(sum(n(S(m2,t-1)+1,:,2))+gamma1+gamma2)-log(2))*ones(1,3)]; 

            end
        end
        s1=[0 Snew(end,t)*ones(1,3) 3-Snew(end,t)];
        s2=[Snew(end,t) 0 Snew(end,t) 3-Snew(end,t) Snew(end,t)];
        for i=1:5
            Saux=S;
            Saux(m1,t)=s1(i);
            Saux(m2,t)=s2(i);
            p(i)=p(i)+mLik_TB(X(:,1:t),Saux(:,1:t),r,chi,tau,nu);
        end
        maximo = max(p);
        p=exp(p-maximo);
        p=p/sum(p);
        pSplit=pSplit+ log(p(find(s1==S(m1,t) & s2==S(m2,t))));

        n(S(m1,t-1)+1,S(m1,t)+1,1)=n(S(m1,t-1)+1,S(m1,t)+1,1)+1;
        n(S(m2,t-1)+1,S(m2,t)+1,2)=n(S(m2,t-1)+1,S(m2,t)+1,2)+1;
    end       
end