cd /export/usuarios01/franrruiz87/simMIMO

T = 2000;
Nt = 6;
Nr = 12;
M = 2;
L = 3;
Niter = 30000;
lHead = 0;
onOffModel = 0;
Nparticles = 3000;
flagParallel = 1;
itCluster = 1;
simId = 3;

tic
simClusterWise(T,Nt,Nr,M,L,Niter,lHead,onOffModel,Nparticles,flagParallel,itCluster,simId)
toc
