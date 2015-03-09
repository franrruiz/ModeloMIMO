function Zest=FBsampling(X,Mvec,H,Nsim,s2g)


p1=0.8;
T=size(X,2);
M=length(Mvec);
Zest=zeros(M,T);
K=max(Mvec);

p0 = [0.8 0.1 0.1];%p0/sum(p0);
ptrans=cell(1,M);
Xhat=cell(1,M);

for m=1:M
    Zest(m,1)=find(mnrnd(1,p0)==1)-1;
    for t=2:T
        Zest(m,t) = find(mnrnd(1,p0)==1)-1;
    end
    Q=3^Mvec(m);
    A=de2bi([0:Q-1]', Mvec(m),3)-1;
    ptrans{m}=zeros(Q);
    
    for q=1:Q
        idx1=bi2de([A(q,2:end),1]+1,3)+1;
        idx0=bi2de([A(q,2:end),0]+1,3)+1;
        idx_1=bi2de([A(q,2:end),-1]+1,3)+1;
       
        if A(q,end)==1 || A(q,end)==-1
            ptrans{m}(q,idx1)=p1/2;
            ptrans{m}(q,idx0)=1-p1;
            ptrans{m}(q,idx_1)=p1/2;
            
        elseif A(q,end)==0
            ptrans{m}(q,idx1)=(1-p1)/2;
            ptrans{m}(q,idx0)=p1;
            ptrans{m}(q,idx_1)=(1-p1)/2;
            
        end
    end 
    ptrans{m}=ptrans{m}./(repmat(sum(ptrans{m},2),1,Q));

end

Zest(find(Zest==2))=-1;

for it=1:Nsim
    for m=1:M
        A=de2bi([0:3^Mvec(m)-1]', Mvec(m),3)-1;
        Xhat{m}=X;
        for m2=[1:m-1 m+1:M]
            Haux=H(:,m2:M:end);
            st=zeros(K,1);
            for t=1:T
                st(1:max(K-t,0))=0;
                st(max(K-t,0)+1:K)=Zest(m2,t+max(K-t,0)-K+1:t);
                Xhat{m}(:,t)= Xhat{m}(:,t)-Haux*st;
            end

        end
        Haux=H(:,m:M:end);
        Haux=Haux(:,K-Mvec(m)+1:end);
        Xopt= Haux*A';
        Q=3^Mvec(m);
        fw= zeros(Q,T);
        bw= zeros(Q,T);

        %%Forward
         t=1;
         fw(:,t)=-1/(2*s2g)*diag((repmat(Xhat{m}(:,t),1,Q)-Xopt)'*(repmat(Xhat{m}(:,t),1,Q)-Xopt))+log(ptrans{m}(find(sum(abs(A),2)==0),:)');
         fw(:,t)=exp(fw(:,t)-max(fw(:,t)));
         fw(:,t)= fw(:,t)/sum( fw(:,t));
         for t=2:T
             fw(:,t)=-1/(2*s2g)*diag((repmat(Xhat{m}(:,t),1,Q)-Xopt)'*(repmat(Xhat{m}(:,t),1,Q)-Xopt))+log(ptrans{m}(:,:)'*fw(:,t-1));
             fw(:,t)=exp(fw(:,t)-max(fw(:,t)));
             fw(:,t)= fw(:,t)/sum(fw(:,t));
         end

         %Backward
         t=T;
         StF=zeros(1,T);
         StF(t)=mnrnd(1,fw(:,t))*[1:Q]';
         for t=T-1:-1:1
            bw(:,t)=fw(:,t).*ptrans{m}(:, StF(t+1)); 
            bw(:,t)= bw(:,t)/sum(bw(:,t));

            StF(t)=mnrnd(1,bw(:,t))*[1:Q]';

         end

         for t=1:T
            sbs=de2bi(StF(t)-1, Mvec(m),3)-1;
            Zest(m,t)=sbs(end);
    %         if t==1 && sbs(end-1)==0
    %             Zest(m,t)=sbs(end);
    %         elseif  sbs(end-1)==Zest(m,t-1)
    %             Zest(m,t)=sbs(end);
    %         else
    %             disp('raro');
    %         end
         end

    end
end