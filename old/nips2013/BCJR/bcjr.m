function Sest=bcjr(X,H,s2g)
T=size(X,2);
p1=0.8;

ptrans=zeros(9,9);
ptrans(1,:)=[p1 0 0 (1-p1)/2 (1-p1)/2 0 0 0 0];
ptrans(2,:)=[p1 0 0 (1-p1)/2 (1-p1)/2 0 0 0 0];
ptrans(3,:)=[p1 0 0 (1-p1)/2 (1-p1)/2 0 0 0 0];
ptrans(4,:)=[0 1-p1 0 0 0 p1/2 0 p1/2 0];
ptrans(5,:)=[0 0 1-p1 0 0 0 p1/2 0 p1/2];
ptrans(6,:)=[0 1-p1 0 0 0 p1/2 0 p1/2 0];
ptrans(7,:)=[0 1-p1 0 0 0 p1/2 0 p1/2 0];
ptrans(8,:)=[0 0 1-p1 0 0 0 p1/2 0 p1/2];
ptrans(9,:)=[0 0 1-p1 0 0 0 p1/2 0 p1/2];
ptrans=kron(ptrans,ptrans);
ptrans=ptrans./(repmat(sum(ptrans,2),1,81));

alph=[0 0; 0 1; 0 -1; 1 0; -1 0; 1 1; 1 -1; -1 1; -1 -1];

A=zeros(81,4);
for i=1:9
    for i2=1:9
        A((i-1)*9+i2,:)=[alph(i,:) alph(i2,:)];
    end
end
Haux=H(:,[3 1 4 2]);
Xopt=Haux*A';

fw=zeros(81,T);
bw=zeros(81,T);
fw(:,1)=ptrans(1,:).*exp(-sum((repmat(X(:,1),1,81)-Xopt).^2,1)/(2*s2g));
fw(:,1)=fw(:,1)/sum(fw(:,1));
bw(:,T)=1;
for t=2:T
    fw(:,t)=(ptrans'*fw(:,t-1)).*exp(-sum((repmat(X(:,t),1,81)-Xopt).^2,1)/(2*s2g))';
    fw(:,t)=fw(:,t)/sum(fw(:,t));
    %aux2 = b_log(:,T-t+2)+log(bw(:,T-t+2));
    %bw(:,T-t+1)= a{l}*(b(:,T-t+2).*bw(:,T-t+2));
    bw(:,T-t+1)= ptrans*(bw(:,T-t+2).*exp(-sum((repmat(X(:,T-t+2),1,81)-Xopt).^2,1)/(2*s2g))');
    bw(:,T-t+1)=bw(:,T-t+1)/sum(bw(:,T-t+1));
end
qt = fw.*bw;
qt=qt./repmat(sum(qt,1),81,1);

qt_red = zeros(9,T);
for t=1:T
    for i=1:9
        indices = find(sum(repmat(alph(i,:),81,1)==A(:,[1 3]),2)==2);
        qt_red(i,t) = sum(qt(indices,t));
    end
end

[val idx]=max(qt_red,[],1);
Sest=zeros(2,T);
for t=1:T
    Sest(:,t)=alph(idx(t), :)';
end

% [val idx]=max(qt,[],1);
% 
% Sest=zeros(2,T);
% for t=1:T
%     Sest(:,t)=A(idx(t), [1 3])';
% end
%errores=sum(sum(abs(simb(:,2:end)-Sest)));