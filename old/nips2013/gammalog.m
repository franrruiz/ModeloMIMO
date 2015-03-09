function y= gammalog(x)
y=gammaln(x);
% y=zeros(size(x));
% for i=1:size(x,1)
%     for j=1:size(x,2)
%         n=x(i,j);
%         if n>171
%             y(i,j)= log(n).*(n-0.5)-n+log(sqrt(2*pi))+log((1+1./(12*n)+1./(288*n.^2)-139./(51840*n.^3)-571./(2488320*n.^4)));
%         else
%             y(i,j)=log(gamma(n));
% 
%         end
%     end
end