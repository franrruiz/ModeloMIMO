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

minRaf = 50; %
idx1=minRaf-1+randi(round(T/2)-minRaf+1);%max(30,randi(T/2));
simb(1,1:idx1)=1-2*(rand(1,idx1)<0.5);
for n=2:N
    idx2 = T+1;
    while(idx2>T)
        idx1 = randi(round(T));
        idx2 = idx1+(minRaf-2)+randi(round(T/2)-minRaf+1);
    end
    simb(n,idx1:idx2)=1-2*(rand(1,idx2-idx1+1)<0.5);
%     if idx2<=T
%         simb(n,idx1:idx2)=1-2*(rand(1,idx2-idx1+1)<0.5);
%     else
%         idx1=T-30-randi(20);
%     	simb(n,idx1:end)=1-2*(rand(1,T-idx1+1)<0.5);
%     end
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
