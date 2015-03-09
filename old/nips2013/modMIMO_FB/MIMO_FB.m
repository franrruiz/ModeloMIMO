function [Zest, rest,Hest,M,nest]= MIMO_FB(X,Zest,nest,rest,chi,tau,nu,alpha,gamma1,gamma2,beta,Nsim)
K=max(rest);
[M T] = size(Zest);
D=size(X,1);
Puntos=[];
Q=3;

Zest(find(Zest==2))=-1;
%% Sample H and Ptrans
[A, Hest, s2x]= sampleAandH(X,Zest,rest,chi,nest,alpha,gamma1,gamma2,tau,nu);


%% Sample Mnew
% sampling the slice variable
minAm=1-max(A(1,1,:));
mu=betarnd(.1,5)*minAm;

varargin=cell(1,3);
varargin{1}=alpha;
varargin{2}=T;
%options=optimset('GradObj','on','Display','none','TolX',1e-6);
options=optimset('Display','none','TolX',1e-12);

flagA=1;
A0no0=[];
Aaux=log(minAm);
while flagA==1
    domain=[-inf,Aaux];
    varargin{3}=exp(Aaux);
    %a=fmincon(@(z)lc_lpdf_grad(z,varargin),log(Aaux)-10,[],[],[],[],[],log(Aaux),[],options);
    a=fminsearchbnd(@(z)lc_lpdf_neg(z,varargin),Aaux-10,-Inf,Aaux,options);
    
    b=Aaux;
    Puntos=Puntos(find(Puntos<b));
    [Aaux Puntos] = ars(@lc_lpdf, a-10, b, domain, 1, Puntos, varargin);
  
    %Aaux=betarnd(alpha,1)*Aaux;
    if Aaux>mu
       A0no0=[A0no0; exp(Aaux)];
    else
       flagA=0;
    end
end

%
Mnew=length(A0no0);
if Mnew>0
    Kold=K;
    rest=[rest 1+poissrnd(beta,1,Mnew)]; 
    K=max(rest);
    Zest=[Zest; zeros(Mnew,T)];
    A=cat(3,A,zeros(Q,Q,Mnew));
    H=zeros(D,K*(M+Mnew));

    for k=1:Kold
        H(:,(K-Kold+k-1)*(M+Mnew)+1:(K-Kold+k-1)*(M+Mnew)+M)= Hest(:,(k-1)*M+1:k*M);    
    end
    
    for m=M+1:M+Mnew
        H(:,(K-rest(m))*(M+Mnew)+m:(M+Mnew):end)=randn(D,rest(m));
    end
    
    Hest=H;

    A(1,1,M+1:end)=1-A0no0;
    A(1,2:end,M+1:end)=reshape(repmat(A0no0,1,Q-1)/(Q-1),1,Q-1,Mnew);
    
    A(2:end,1,M+1:end)=reshape(betarnd(gamma1,gamma2,Q-1,Mnew),Q-1,1,Mnew);
    A(2:end,2:end,M+1:end)=repmat(1-A(2:end,1,M+1:end),[1,Q-1,1])/(Q-1);
    
    M=M+Mnew;
    
end

%% FB
Zest=FBsampling(X,rest,Hest,A,Zest,Nsim,s2x);

mNOempt=find(sum(Zest~=0,2)>0);
Zest=Zest(mNOempt,:);
A=A(:,:,mNOempt);
rest=rest(mNOempt);


H=zeros(D,K*length(mNOempt));
for k=1:K
    H(:,(k-1)*length(mNOempt)+1:k*length(mNOempt))= Hest(:,(k-1)*M+mNOempt);
end
Knew=max(rest);
M=length(mNOempt);
Hest=H(:,(K-Knew)*M+1:end);
K=Knew;


if isempty(Zest)
    Zest=zeros(1,T);
    M=1;
    rest=poissrnd(beta)+1;
    Hest=zeros(D,rest);
    A=zeros(Q,Q,1);
end

% Computing nest
Zest(find(Zest==-1))=2;
nest=zeros(Q,Q,M);
for m=1:M
    for q=1:Q
        for k=1:Q
            nest(q,k,m)=sum([0 Zest(m,1:end-1)]==q-1 & Zest(m,:)==k-1);
        end
    end
end