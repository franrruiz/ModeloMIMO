%% Modelo: Canal MIMO variante en el Tiempo
% T=200;
% %%Definimos el canal
% N=1;        % entradas
% L = 20;		% salidas
% s2g=1;
% M=1;
% s2h=10;

H=sqrt(s2h)*randn(L,N*M);
simb=zeros(N,T);

idx1=randint(1,1,[30 round(T/2)]);
simb(1,1:idx1)=1-2*(rand(1,idx1)<0.5);
for n=2:N
   idx1=round(T/(2*N))*(n-1)+randint(1,1,[1 round(T/(6*N))]);
   idx2= idx1+randint(1,1,[30 T]);
   if idx2<=T
      simb(n,idx1:idx2)=1-2*(rand(1,idx2-idx1+1)<0.5);
   else
       num=randint(1,1,[30 60]);
       simb(n,T-num+1:T)=1-2*(rand(1,num)<0.5);
   end
end
simb=[zeros(N,M-1) simb];

X=zeros(L,T);

for t=1:T
   y=simb(:,t:t+M-1); 
   X(:,t)=H*y(:)+sqrt(s2g)*randn(L,1);
end
