N=3;
T=400;
M=3;
simb=zeros(N,T);

minRaf = 80; %
idx1=minRaf-1+randint(1,1,[1 round(T/2)-minRaf+1]);%max(30,randi(T/2));
simb(1,1:idx1)=1-2*(rand(1,idx1)<0.5);
for n=2:N
    idx2 = T+1;
    while(idx2>T)
        idx1 = randint(1,1,[1 round(T)]);
        idx2 = idx1+(minRaf-2)+randint(1,1,[1 round(T/2)-minRaf+1]);
    end
    simb(n,idx1:idx2)=1-2*(rand(1,idx2-idx1+1)<0.5);
    
end
simb=[zeros(N,M-1) simb];

save('symbols.mat','simb');