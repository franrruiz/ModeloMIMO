%% Modelo: Canal MIMO variante en el Tiempo
% T=200;
% %%Definimos el canal
% N=1;        % entradas
% L = 20;		% salidas
% s2g=1;
% M=1;
% s2h=10;

M=max(Mvec);
H=sqrt(s2h)*randn(L,N*M);
simb=zeros(N,T);
if Mvec(1)<M
      H(:,1:N:1+N*(M-Mvec(1))-1)=0; 
 end

idx1=max(30,randint(1,1,[1 round(T/2)]));
simb(1,1:idx1)=1-2*(rand(1,idx1)<0.5);
for n=2:N
   idx1=randint(1,1,[1 round(T)]);
   idx2= idx1+29+randint(1,1,[1 round(T/2-29)]);
   if idx2<=T
    simb(n,idx1:idx2)=1-2*(rand(1,idx2-idx1+1)<0.5);
   else
    idx1=T-30-randint(1,1,[1 20]);
    simb(n,idx1:end)=1-2*(rand(1,T-idx1+1)<0.5);
   end
   if Mvec(n)<M
      H(:,n:N:n+N*(M-Mvec(n))-1)=0; 
   end
end
simb=[zeros(N,M-1) simb];

X=zeros(L,T);

for t=1:T
   y=simb(:,t:t+M-1); 
   X(:,t)=H*y(:)+sqrt(s2g)*randn(L,1);
end
