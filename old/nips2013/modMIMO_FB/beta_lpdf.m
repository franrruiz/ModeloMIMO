function pdfz=beta_lpdf(z,varargin)

alpha=varargin{1}{1};
beta=varargin{1}{2};

if alpha==1 && beta~=1
    pdfz=(beta-1)*log(1-z);
elseif beta==1 && alpha~=1
    pdfz=(alpha-1)*log(z);
elseif alpha==1 && beta==1
    pdfz=zeros(size(z));
else
    pdfz=(alpha-1)*log(z)+(beta-1)*log(1-z);
end

end