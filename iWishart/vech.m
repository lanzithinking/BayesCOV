% Gauss vech function
% vectorize a lower triangular matrix in the chosen order
% -------------------------------------------------------
% Shiwei Lan @ CalTech, 2017
% __________________________

function v=vech(L,order)
if ~exist('order','var')
    order='col';
end

D=size(L,1);

if any(strfind('column',order))
    ind=tril(true(D));
    v=L(ind);
elseif any(strfind('row',order))
    ind=triu(true(D));
    R=L';
    v=R(ind);
else
    warning('Wrong order!');
    v=nan(D*(D+1)/2,1);
end

end