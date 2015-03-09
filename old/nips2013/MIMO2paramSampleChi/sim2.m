randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

T=400;
%%Definimos el canalc
L = 10;	% salidas
SNR=6;
s2g= 10^(-SNR/10);
Mvec=[3 2 1];
N=length(Mvec);

s2h=1;

generaMIMOrafagas2;
chi=1;%1/s2g;%2*s2g;
tau=1;
nu=.1;
alpha1=.1;
alpha2=1;
gamma1=.1;
gamma2=10;
beta=1;
pii=0.2;
lambda=1;
D=L;


Nsim=100;



[Zest,Hest,MestIt,pX_ZIt,chiIt,rest,nest] = MIMO_aprender_H_s2x(X,chi,tau,nu,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, 1);

% save('Sim1Final.mat');
% 
% exit
