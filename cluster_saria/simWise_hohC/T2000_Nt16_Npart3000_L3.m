cd /export/usuarios01/franrruiz87/simMIMO

T = 2000;
Nt = 16;
Nr = 12;
M = 2;
L = 3;
Niter = 50000;
lHead = 0;
onOffModel = 0;
Nparticles = 3000;
blockNtSize = 9;
stepDB_thr = 40000;
stepDB_ini = 0.001;
stepDB_end = 0.008;
flagParallel = 1;
itCluster = 1;
simId = 3;

%tic
simClusterWise(T,Nt,Nr,M,L,Niter,lHead,onOffModel,Nparticles,blockNtSize,stepDB_thr,stepDB_ini,stepDB_end,flagParallel,itCluster,simId);
%toc
