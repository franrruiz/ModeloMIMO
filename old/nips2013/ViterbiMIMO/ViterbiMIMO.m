%function ViterbiMIMO(X,H, Mvec);
T=size(X,2);
Nst=3^(sum(Mvec-1));
N=length(Mvec);


Tsimb=de2bi([0:Nst-1]',sum(Mvec-1),3);
Tsimb=Tsimb-1;

Nss=3^(sum(Mvec>1));
Tss=zeros(Nst,Nss);

for s=1:Nst
    inter=1:Nst;
    
    for tx=1:(sum(Mvec>1))
        un=union([],find(sum(repmat([-1 Tsimb(s,sum(Mvec(1:tx-1))-tx+2:sum(Mvec(1:tx))-tx-1)],Nst,1)==Tsimb(:,sum(Mvec(1:tx-1))-tx+2:sum(Mvec(1:tx))-tx),2)==Mvec(tx)-1));
        un=union(un,find(sum(repmat([0 Tsimb(s,sum(Mvec(1:tx-1))-tx+2:sum(Mvec(1:tx))-tx-1)],Nst,1)==Tsimb(:,sum(Mvec(1:tx-1))-tx+2:sum(Mvec(1:tx))-tx),2)==Mvec(tx)-1));
        un=union(un,find(sum(repmat([1 Tsimb(s,sum(Mvec(1:tx-1))-tx+2:sum(Mvec(1:tx))-tx-1)],Nst,1)==Tsimb(:,sum(Mvec(1:tx-1))-tx+2:sum(Mvec(1:tx))-tx),2)==Mvec(tx)-1));
        inter=intersect(inter,un);
    end
    
    Tss(s,:)=inter;
end
        



Tm=zeros(Nst,T);
Ts=zeros(Nst,N,T);

Np=3^(sum(Mvec==1));
Tnp=de2bi([0:Np-1]',sum(Mvec==1),3);
Tnp=Tnp-1;

s=find(sum(abs(Tsimb),2)==0);
t=1;
Kmax=max(Mvec);

for ss=1:Nss
    met=inf;
    S=zeros(N*Kmax,1);
    for tx=1:(sum(Mvec>1))
        S(N*(Kmax-1)+tx)=Tsimb(Tss(s,ss),sum(Mvec(1:tx-1))-tx+2);
    end
    
    for p=1:Np
        S(end-sum(Mvec==1)+1:end)=Tnp(p,:);
       
        Lk=sum((X(:,t)-H*S).^2);
        if Lk<met
            met=Lk;
            Sopt=S;
        end
    end
    Ts(Tss(s,ss),:,t)=Sopt(end-N+1:end);
    Tm(Tss(s,ss),t)=met;
end


for t=1:T
    for s=1:Nst
        for ss=1:Nss
            met=1000;
            for p=1:Np
                
                
                
                
            end
        end
    end
end