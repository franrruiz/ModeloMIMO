cd /export/usuarios01/franrruiz87/simMIMO

T = 1000;
Nt = 10;
Nr = 20;
M = 2;
Ltrue = 1;
L = 1;
SNR = -3;
Niter = 10000;
lHead = 0;
onOffModel = 0;
Nparticles = 300;
flagParallel = 1;
itCluster = 1;
simId = 21;

tic
simClusterSynGenie(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,Nparticles,flagParallel,itCluster,simId);
toc
