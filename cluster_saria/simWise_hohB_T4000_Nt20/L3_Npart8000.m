cd /export/usuarios01/franrruiz87/simMIMO

T = 4000;
Nt = 20;
Nr = 11;
M = 2;
Ltrue = 3;
L = 5;
SNR = 10;
Niter = 10000;
lHead = 0;
onOffModel = 0;
Nparticles = 8000;
flagParallel = 1;
itCluster = 1;
simId = 2;

tic
simClusterWise(T,Nt,Nr,M,Ltrue,L,SNR,Niter,lHead,onOffModel,Nparticles,flagParallel,itCluster,simId);
toc
