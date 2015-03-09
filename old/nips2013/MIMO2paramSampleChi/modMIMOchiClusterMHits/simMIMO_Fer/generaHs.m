%% Modelo: Canal MIMO variante en el Tiempo
% T=200;
% %%Definimos el canal
% N=1;        % entradas
% L = 20;		% salidas
% s2g=1;
% M=1;
% s2h=10;
L=10;
Mvec=[3 2 1];
M=max(Mvec);
s2h=1;
N=length(Mvec);
T=400;

for it=1:20
    H=sqrt(s2h)*randn(L,N*M);

    if Mvec(1)<M
          H(:,1:N:1+N*(M-Mvec(1))-1)=0; 
    end

    for n=2:N
        if Mvec(n)<M
            H(:,n:N:n+N*(M-Mvec(n))-1)=0; 
        end
    end
    save(['H_it' num2str(it) '.mat'],'H');
end

