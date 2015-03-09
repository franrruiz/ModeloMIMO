function [Snew nnew pMerge flagSigno]=mergeChains(S,n,m1,m2,t1,t2)
T=size(S,2);
Snew=S([1:min(m1,m2)-1 min(m1,m2)+1:max(m1,m2)-1 max(m1,m2)+1:end],:);
flagSigno=0;
if rand>0.5
    idx=find(S(m2,:)>0);
    S(m2,idx)=3-S(m2,idx);
    flagSigno=1;
end

Snew=[Snew; zeros(1,T)];
idxI=find(S(m1,:)==S(m2,:));
idx0=find((S(m1,:)==0 & S(m2,:)~=0) | (S(m1,:)~=0 & S(m2,:)==0));
idxD=find(S(m1,:)~=0 & S(m2,:)~=0 & S(m1,:)~=S(m2,:));
idxD=setdiff(idxD,[t1 t2]);

Snew(end,idxI)= S(m1,idxI);
Snew(end,idx0)= S(m1,idx0)+S(m2,idx0);
Snew(end,idxD)=(rand(1,length(idxD))>0.5)+1;
Snew(end,t1)=S(m1,t1);
Snew(end,t2)=S(m2,t2);

nnew=n(:,:,[1:min(m1,m2)-1 min(m1,m2)+1:max(m1,m2)-1 max(m1,m2)+1:end]);
nnew =cat(3,nnew,zeros(3,3));

m=size(Snew,1);
for q=1:3
    for k=1:3
        nnew(q,k,m)=sum([0 Snew(m,1:end-1)]==q-1 & Snew(m,:)==k-1);
    end
end

pMerge=-length(idxD)*log(2)-log(2);
