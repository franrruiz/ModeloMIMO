function [Sest SeqEst nest flagAccept] = ep_main(data,samples,hyper,param)

if(param.L~=1)
    error('EP not implemented for L~=1');
end

Nt = size(samples.Z,1);
samplesAux = samples;
samplesAux.Z = zeros(size(samples.Z));
samplesAux.seq = zeros(size(samples.seq));

if(param.flag0)
    constellation = [0 param.constellation];
else
    constellation = param.constellation;
end

H = [real(samples.H) -imag(samples.H);
     imag(samples.H) real(samples.H)];

Es = sqrt(mean(abs(constellation).^2));

log_qZ = 0;
log_qZnew = 0;
for t=1:param.T
    % Run EP:
    yt = [real(data.obs(:,t)); imag(data.obs(:,t))];
    [mm,SS] = estimatorEP(constellation,H,yt,Es,samples.s2y,param.ep.eps,param.ep.beta,param.ep.Niter);
    
    % Obtain the component-wise Gaussian means and variances:
    aux = mm;
    ss = diag(SS);
    ss = [ss(1:Nt) ss(Nt+1:2*Nt)];
    
    % Evaluate the Gaussian in the points of the constellation:
    distancias = zeros(Nt,length(constellation));   % Size = [Nt x Q]
    for jj=1:length(constellation)
        distancias(:,jj) = -sum(((aux-repmat([real(constellation(jj)) imag(constellation(jj))],Nt,1)).^2)./(2*ss),2);
    end
    distancias = exp(distancias-repmat(max(distancias,[],2),1,length(constellation)));
    distancias = distancias./repmat(sum(distancias,2),1,length(constellation));

    % Sample constellation points from the normalized distribution:
    idx = randmult2(distancias.').';   % Size = [Nt x 1]
    
    % Compute q(Z) and q(Z*) for instant t:
    idx_qZnew = sub2ind(size(distancias),(1:Nt).',idx);
    log_qZnew = log_qZnew+sum(log(distancias(idx_qZnew)));

    idx_qZ = sub2ind(size(distancias),(1:Nt).',(samples.seq(:,t)+param.flag0));
    log_qZ = log_qZ+sum(log(distancias(idx_qZ)));
    
    % Build samplesAux.Z and samplesAux.seq:
    if(param.flag0)
        idxNZ = (idx>1);
        samplesAux.Z(:,t) = constellation(idx);
        samplesAux.seq(idxNZ,t) = idx(idxNZ)-1;
    else
        samplesAux.Z(:,t) = constellation(idx);
        samplesAux.seq(:,t) = idx;
    end
end

% Create struct samplesAux:
samplesAux.nest = zeros(2,2,size(samplesAux.Z,1));
for m=1:size(samplesAux.Z,1)
    % From 0 to 0
    samplesAux.nest(1,1,m) = sum([0 samplesAux.seq(m,1:end-1)]==0 & samplesAux.seq(m,:)==0);
    % From 0 to active
    samplesAux.nest(1,2,m) = sum([0 samplesAux.seq(m,1:end-1)]==0 & samplesAux.seq(m,:)~=0);
    % From active to 0
    samplesAux.nest(2,1,m) = sum([0 samplesAux.seq(m,1:end-1)]~=0 & samplesAux.seq(m,:)==0);
    % From active to active
    samplesAux.nest(2,2,m) = sum([0 samplesAux.seq(m,1:end-1)]~=0 & samplesAux.seq(m,:)~=0);
end

% Compute acceptance probability:
pA = compute_llh(data,samplesAux,hyper,param)-compute_llh(data,samples,hyper,param) ...
     +compute_pZ(data,samplesAux,hyper,param)-compute_pZ(data,samples,hyper,param) ...
     +log_qZ-log_qZnew;
if(rand(1)<exp(pA))
    % Accept
    flagAccept = 1;
    Sest = samplesAux.Z;
    SeqEst = samplesAux.seq;
    nest = samplesAux.nest;
else
    % Reject
    flagAccept = 0;
    Sest = samples.Z;
    SeqEst = samples.seq;
    nest = samples.nest;
end
