function pdfz=lc_lpdf_neg(z,varargin)

alpha=varargin{1}{1};
T=varargin{1}{2};
cant=varargin{1}{3};

tvec=(1:T)';

pdfz=-(alpha*z+alpha*sum(repmat(1-exp(z),T,1).^repmat(tvec,1,length(z))./repmat(tvec,1,length(z)),1)+T*log(1-exp(z)));

end