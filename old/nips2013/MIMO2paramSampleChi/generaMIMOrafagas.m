%% Modelo: Canal MIMO variante en el Tiempo
% T=200;
% %%Definimos el canal
% N=1;        % entradas
% L = 20;		% salidas
% s2g=1;
% M=1;
% s2h=10;

H=sqrt(s2h)*randn(L,N*M);

simb=1-2*(rand(N,T)<0.5);
simb=[zeros(N,M-1) simb];
simb(1,1:randint(1,1,[1 T/2]))=0;
for n=2:N
   if n<N
    simb(n,1:randint(1,1,[1 T/2]))=0;
   else
   t1=randint(1,1,[1 T]);
   dur=randint(1,1,[1 T/2]);
   if t1+dur>T+M-1
       simb(n,t1:end)=0;
   else
       simb(n, t1:t1+dur)=0;
   end
   end
end
X=zeros(L,T);

for t=1:T
   y=simb(:,t:t+M-1); 
   X(:,t)=H*y(:)+sqrt(s2g)*randn(L,1);
end
