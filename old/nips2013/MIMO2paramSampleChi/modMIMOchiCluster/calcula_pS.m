function pS= calcula_pS(n,alpha1,alpha2,gamma1,gamma2,T)
M=size(n,3);

pS=M*log(alpha1*alpha2)-alpha2*alpha1*sum(1./((1:T)+alpha1-1));

for m=1:M
    pS=pS-sum(n(1,2:end,m))*log(2)+gammaln(alpha1+n(1,1,m))+gammaln(sum(n(1,2:end,m)))-gammaln(alpha1+sum(n(1,:,m))) + 2*gammaln(gamma1+gamma2)-2*gammaln(gamma1)-2*gammaln(gamma2)-...
         sum(n(2,2:end,m))*log(2)+gammaln(gamma1+n(2,1,m))+gammaln(gamma2+sum(n(2,2:end,m)))-gammaln(gamma1+gamma2+sum(n(2,:,m)))-...
         sum(n(3,2:end,m))*log(2)+gammaln(gamma1+n(3,1,m))+gammaln(gamma2+sum(n(3,2:end,m)))-gammaln(gamma1+gamma2+sum(n(3,:,m)));
end
end