function pX_Zlog = mLik_TB(X,S,r,s2x,s2h) %%TB=tapas bar
%% Likelihood computation
%% Input
% X is the discrete observation matrix, size D (dimensionality)xT(number of observatios)
% Z IBP matrix
% s2x variance of the observations
% s2phi variance of the same state among the chains
% s20 variance of 
%% Output
% pX_Zlog likelihood (p(X|Z), where Phi parameter have been integrated out)

[D T]=size(X);
M=size(S,1);
K=max(r);
Q=3;
S=S';
S(find(S==2))=-1;
Z=zeros(T,M*K);
R=zeros(M*K,M*K);

for k=1:K
    Z(k:end,(k-1)*M+1:k*M)=S(1:T-k+1,:);
    R((k-1)*M+1:k*M,(k-1)*M+1:k*M)=diag(r>=k);
end
ZR=Z*R;
W=ZR'*ZR+s2x/s2h*eye(M*K);
Sigma=eye(T)-ZR*(W\ZR') ;


pX_Zlog=-T*D/2*log(2*pi)-T*D/2*log(s2x)+D/2*log(det(Sigma))-1/(2*s2x)*trace(X*Sigma*X');
