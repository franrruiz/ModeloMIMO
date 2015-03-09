randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

T=150;
%%Definimos el canal
N=4;        % entradas
L = 30;	% salidas
s2g=.1;
Mvec=[5 1 3 2];

s2h=1;

generaMIMOrafagas2;
chi=2*s2g;
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


Nsim=200;



[Zest,Hest,MestIt,pX_ZIt,rest,nest] = MIMO_aprender_H_s2x(X,chi,tau,nu,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, 1);

% save('Sim1Final.mat');
% 
% exit
