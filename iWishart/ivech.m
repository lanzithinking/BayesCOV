% inverse Gauss vech function
% restore the lower triangular matrix from a vector with the chosen order
% -------------------------------------------------------
% Shiwei Lan @ CalTech, 2017
% __________________________

function L=ivech(v,order)
if ~exist('order','var')
    order='col';
end

l=length(v);
D=ceil((sqrt(1+8*l)-1)/2);

L=zeros(D);
if any(strfind('column',order))
    ind=tril(true(D));
    L(ind)=v;
elseif any(strfind('row',order))
    ind=triu(true(D));
    L(ind)=v;
    L=L';
else
    warning('Wrong order!');
end

end