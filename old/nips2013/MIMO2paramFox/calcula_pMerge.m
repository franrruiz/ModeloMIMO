function pMerge=calcula_pMerge(Snew,m1,t1,t2)

m2=size(Snew,1);

idxD=find(Snew(m1,:)~=0 & Snew(m2,:)~=0 & Snew(m1,:)~=Snew(m2,:));
idxD=setdiff(idxD,[t1 t2]);

pMerge=-length(idxD)*log(2)-log(2);