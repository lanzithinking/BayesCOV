% inverse Gauss vech function
% restore the lower triangular matrix from a vector with the chosen order
% extended to time series data
% -------------------------------------------------------
% Shiwei Lan @ CalTech, 2017
% __________________________

function L=ivechx(v,order)
if ~exist('order','var')
    order='col';
end

l=size(v,2); % l : (N,D*(D+1)/2)
D=ceil((sqrt(1+8*l)-1)/2);

L=zeros(size(v,1),D,D);
if any(strfind('column',order))
    ind=tril(true(D));
    L(:,ind)=v;
elseif any(strfind('row',order))
    ind=triu(true(D));
    L(:,ind)=v;
    L=permute(L,[1,3,2]);
else
    warning('Wrong order!');
    L=nan(size(v,1),D,D);
end

end