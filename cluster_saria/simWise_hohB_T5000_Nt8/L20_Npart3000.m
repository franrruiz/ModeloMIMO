cd /export/usuarios01/franrruiz87/simMIMO

T = 5000;
Nt = 8;
Nr = 11;
M = 2;
Ltrue = 20;
L = 20;
SNR = 10;
Niter = 10000;
lHead = 0;
onOffModel = 0;
Nparticles = 3000;
flagParallel = 1;
itCluster = 1;
simId = 2;

tic
simClusterWise(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,Nparticles,flagParallel,itCluster,simId);
toc
