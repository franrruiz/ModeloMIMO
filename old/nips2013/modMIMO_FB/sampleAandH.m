function [A, H, s2x]= sampleAandH(X,Zest,rest,chi,nest,alpha,gamma1,gamma2,tau,nu)
[M T]=size(Zest);
D=size(X,1);
Q=3;
%% Sampling A
A=zeros(Q,Q,M);
if M>1
    % Sample a(1)
    m=1;
    a0no0 = betarnd(alpha+sum(nest(1,2:end,m)), nest(1,1,m)+1);
    A(1,1,m)= 1-a0no0;
    A(1,2:end,m)=a0no0/(Q-1);
    for q=2:Q
        A(q,1,m)=betarnd(gamma1+nest(q,1,m), gamma2+sum(nest(q,2:end,m)));
        A(q,2:end,m)=(1-A(q,1,m))/(Q-1);
    end
    % Sample a(2:M-1)
    for m=2:M
        %a0no0=betaTruncrnd(alpha+sum(nest(1,2:end,m)), nest(1,1,m)+1,1-A(1,1,m+1), 1-A(1,1,m-1));
        varargin=cell(1,2);
        varargin{1}=alpha+sum(nest(1,2:end,m));
        varargin{2}=nest(1,1,m)+1;
        domain = [0 1-A(1,1,m-1)];
        [a0no0 Puntos] = ars(@beta_lpdf, domain(1), domain(2), domain, 1, [], varargin);
        A(1,1,m)= 1-a0no0;
        A(1,2:end,m)=a0no0/(Q-1);
        for q=2:Q
            A(q,1,m)=betarnd(gamma1+nest(q,1,m), gamma2+sum(nest(q,2:end,m)));
            A(q,2:end,m)=(1-A(q,1,m))/(Q-1);
        end
    end
else   % M==1
    % Sample a(1)
    m=1;
    a0no0 = betarnd(alpha+sum(nest(1,2:end,m)), nest(1,1,m)+1);
    A(1,1,m)= 1-a0no0;
    A(1,2:end,m)=a0no0/(Q-1);
    for q=2:Q
        A(q,1,m)=betarnd(gamma1+nest(q,1,m), gamma2+sum(nest(q,2:end,m)));
        A(q,2:end,m)=(1-A(q,1,m))/(Q-1);
    end
end
    

%% Sample s2x
K=max(rest);
S=Zest';
%S(find(S==2))=-1;
Z=zeros(T,M*K);
R=zeros(M*K,M*K);
for k=1:K
    Z(k:end,(k-1)*M+1:k*M)=S(1:T-k+1,:);
    R((k-1)*M+1:k*M,(k-1)*M+1:k*M)=diag(rest>=k);
end
ZR=Z*R;
W_1=ZR'*ZR+chi*eye(M*K);
Sigma_1=eye(T)-ZR*(W_1\(ZR'));

s2x=1/(gamrnd(T*D/2+tau,1/(nu+0.5*trace(X*Sigma_1*X'))));


%% Sampling H

SigmaH_1=W_1/s2x;
MH=ZR'*X';

% % Check that Sigma_P is symmetric and positive semidefinite
% if (sum(eig(SigmaP_1)<=0)>0)
%     SigmaP_1 = SigmaP_1+eps*eye(size(SigmaP_1));
%     warning('MATLAB:sampleAandPhi:NonPositiveSemidefiniteMatrix','Sigma_P is not positive semidefinite, some constant has been added');
% end

Haux=mvnrnd2((SigmaH_1\MH)',inv(SigmaH_1),D);
H = zeros(D,M*K);
for k=1:K
    H(:,(K-k)*M+1:(K-k+1)*M) = Haux(:,(k-1)*M+1:k*M); 
end

for m=1:M
    H(:,m:M:(K-rest(m))*M)=0;
end