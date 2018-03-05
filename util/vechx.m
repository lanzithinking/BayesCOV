% Gauss vech function
% vectorize a lower triangular matrix in the chosen order
% extended to time series data
% -------------------------------------------------------
% Shiwei Lan @ CalTech, 2017
% __________________________

function v=vechx(L,order)
if ~exist('order','var')
    order='col';
end

D=size(L,2); % L : (N,D,D)

if any(strfind('column',order))
    ind=tril(true(D));
    v=L(:,ind);
elseif any(strfind('row',order))
    ind=triu(true(D));
    R=permute(L,[1,3,2]);
    v=R(:,ind);
else
    warning('Wrong order!');
    v=nan(size(L,1),D*(D+1)/2);
end

end