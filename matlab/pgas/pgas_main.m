function [Sest SeqEst] = pgas_main(data,samples,hyper,param)

%% Extract parameters from structs
Nt = size(samples.H,2);
M = param.pgas.Niter;
N_PF = param.pgas.N_PF;
N_PG = param.pgas.N_PG;
returnN = param.pgas.returnNsamples;

Q = length(param.constellation);
if(param.flag0)
    constellation = [0 param.constellation];
else
    constellation = [param.constellation];
end

%% Call the pgas function
if(isempty(samples.seq))
    flagPG = 0;
else
    flagPG = 1;
end
X_PG = pgas_C(Nt,param.Nr,N_PF,N_PG,M,param.T,param.L,length(constellation),constellation,data.obs,flagPG,samples.seq,samples.H,samples.s2y,samples.am,samples.bm,param.header,length(param.header),param.pgas.particles,param.onOffModel);

%% Return last obtained samples in a [Nt x T x M] matrix
SeqEst = permute(X_PG(:,M-returnN+1:M,:),[1 3 2]);
idxN0 = find(SeqEst~=0);
Sest = zeros(Nt,param.T,returnN);
Sest(idxN0) = constellation(SeqEst(idxN0)+param.flag0);
