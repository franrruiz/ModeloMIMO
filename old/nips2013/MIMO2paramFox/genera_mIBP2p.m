% Genera MIMO dos param.

T=100;
alpha1=.1;
alpha2=1;
gamma1=.1;
gamma2=10;


%t=1
M=Poissrnd(alpha2);
Z=zeros(M,T);
Z(:,1)=(rand(M,1)>0.5)+1;
Q=3;
n=zeros(Q,Q,M);
for m=1:M
    n(1,Z(m,1)+1,m)=1;
end



for t=2:T
    for m=1:M
        if Z(m,t-1)==0
            if rand<sum(n(1,2:end,m))/(alpha1+sum(n(1,:,m)))
                Z(m,t)=(rand>0.5)+1;
                n(1,Z(m,t)+1,m)=n(1,Z(m,t)+1,m)+1;
            else
                n(1,Z(m,t)+1,m)=n(1,Z(m,t)+1,m)+1;
            end
        else 
            if rand<(gamma2+sum(n(Z(m,t-1)+1,2:end,m)))/(gamma1+gamma2+sum(n(Z(m,t-1)+1,:,m)))
                Z(m,t)=(rand>0.5)+1;
                n(Z(m,t-1)+1,Z(m,t)+1,m)=n(Z(m,t-1)+1,Z(m,t)+1,m)+1;
            else
                n(Z(m,t-1)+1,Z(m,t)+1,m)=n(Z(m,t-1)+1,Z(m,t)+1,m)+1;
            end
        end
    end
    Mnew=Poissrnd(alpha1*alpha2/(alpha1+t-1));
    Z=[Z; zeros(Mnew,T)];
    n=cat(3,n,zeros(Q,Q,Mnew));
    for m=M+1:M+Mnew
        Z(m,t)=(rand>0.5)+1;
        n(1,Z(m,t)+1,m)=1;
        n(1,1,m)=t-1;
    end
    M=M+Mnew;
end

                
                
                
                
    