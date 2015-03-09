randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

T=150;
%%Definimos el canal
N=3;        % entradas
L = 20;	% salidas
s2g=1;
M=1;

s2h=1;

generaMIMOrafagas;
s2x=2*s2g;
alpha1=.1;
alpha2=1;
gamma1=.1;
gamma2=10;
beta=1;
pii=0.2;
lambda=1;
D=L;


Nsim=100;



[Zest,Hest,MestIt,pX_ZIt,rest,nest] = MIMO_aprender(X,s2x,s2h,alpha1, alpha2,gamma1,gamma2,beta,pii,lambda,Nsim, 1);

% save('Sim1Final.mat');
% 
% exit
