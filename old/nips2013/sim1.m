randn('seed',round(sum(1e5*clock)));
rand('seed',round(sum(1e5*clock)));

T=100;
%%Definimos el canal
N=3;        % entradas
L = 20;	% salidas
s2g=5;
M=1;

s2h=1;

generaMIMOrafagas;
s2x=s2g;
alpha=1;
gamma1=.1;
gamma2=10;
beta=1;
pii=0.2;
lambda=1;
D=L;


Nsim=100;



[Zest,Hest,MestIt,pX_ZIt,rest,nest] = MIMO_aprender(X,s2x,s2h,alpha,gamma1,gamma2,beta,pii,lambda,Nsim, 1);

% save('Sim1Final.mat');
% 
% exit
