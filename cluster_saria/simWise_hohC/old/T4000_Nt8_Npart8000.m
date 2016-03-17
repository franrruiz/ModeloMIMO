cd /export/usuarios01/franrruiz87/simMIMO

T = 4000;
Nt = 8;
Nr = 12;
M = 2;
L = 5;
Niter = 25000;
lHead = 0;
onOffModel = 0;
Nparticles = 8000;
flagParallel = 1;
itCluster = 1;
simId = 3;

tic
simClusterWise(T,Nt,Nr,M,L,Niter,lHead,onOffModel,Nparticles,flagParallel,itCluster,simId)
toc
