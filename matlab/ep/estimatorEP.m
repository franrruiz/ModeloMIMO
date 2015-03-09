function [mm,Sigma]=estimatorEP(constellation,H,ymode,Es,varruido,eps,Beta,NN)
%Help
%constellation the M points 
%H channel Matrix
%ymode observe vector
%Es energy of constellation
%varruido variance of the noise
%Eps parameter of minimun variance
%Beta  parameter for steps
%NN number of iterations


n=size(H,2);
sigmaruido=varruido;
b=(sigmaruido^-1)*transpose(H)*ymode; 
A=(sigmaruido^-1)*transpose(H)*H;  

ntx=n;
gamma=zeros(ntx,1); % init algorithm
Lambda=(sigmaruido/Es)*ones(ntx,1);



tmpo=real(constellation);re=sort(tmpo(:));re(re((1:end-1)')==re((2:end)'))=[];
tmpo=imag(constellation);im=sort(tmpo(:));im(im((1:end-1)')==im((2:end)'))=[];

if length(im)==1 %check 2QAM
    im=im*zeros(size(re));
end
vme=[(re*ones(1,n/2))';(im*ones(1,n/2))']; 


Sigma= (A+diag(Lambda))\eye(size(A));
sigma=diag(Sigma); % (12)
media=Sigma*(b+gamma);% (11)


h_m=zeros(ntx,NN);
t_m=zeros(ntx,NN);
mediap=zeros(ntx,NN);
sigmap=zeros(ntx,NN);
for it=1:NN 
    
    
    
    h_m(:,it)=max((sigma(:,it).^-1 -Lambda(:,it)).^-1,eps); % 14
    t_m(:,it)=(h_m(:,it)).*((media(:,it).*sigma(:,it).^-1)-gamma(:,it));% 15.
    
    
    
    tmp=normpdf(vme,t_m(:,it)*ones(1,length(vme(1,:))),sqrt(h_m(:,it))*ones(1,length(vme(1,:))));
    Z_i=max(sum(tmp,2),1e-20); % 16 y 25,
    tmpn=tmp./(Z_i*ones(1,length(vme(1,:)))); %  1
    
    
    mediap(:,it)=sum((vme).*(tmpn),2);
    
   
    sigmap(:,it)=max(sum((tmpn).*((vme-mediap(:,it)*ones(1,length(vme(1,:)))).^2),2),eps);
    
   
    % Check for wrong movement
    Lambda(:,it+1)=((1./sigmap(:,it))-(1./h_m(:,it))).*((1./sigmap(:,it))-(1./h_m(:,it))>0)+...
    Lambda(:,it).*((1./sigmap(:,it))-(1./h_m(:,it))<0);
    gamma(:,it+1)=(mediap(:,it)./sigmap(:,it))- (t_m(:,it)./h_m(:,it));
    gamma(:,it+1)=gamma(:,it+1).*((1./sigmap(:,it))-(1./h_m(:,it))>0)+((1./sigmap(:,it))-(1./h_m(:,it))<0).*gamma(:,it);
    
    %OPTION FILTER
    Lambda(:,it+1)=Beta*Lambda(:,it+1)+(1-Beta)*Lambda(:,it);
    gamma(:,it+1)=Beta*gamma(:,it+1)+(1-Beta)*gamma(:,it);
     
    

    Sigma=(A+diag(Lambda(:,it+1)))\eye(size(A));
    sigma(:,it+1)=diag(Sigma);
    media(:,it+1)=Sigma*(b+gamma(:,it+1));

end

mm=reshape(media(:,end),n/2,[]);
%mm2=(mm(:,1)+1j*mm(:,2));
