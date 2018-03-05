% Linear index in vech from multiple subscripts 
% ---------------------------------------------
% Shiwei Lan @ CalTech, 2017
% __________________________

function IND = sub2ind_vech(SQDIM,I,J,order)
if ~exist('order','var')
    order='col';
end

idx_mat = ivech(1:SQDIM*(SQDIM+1)/2,order);
IND = idx_mat(sub2ind([SQDIM,SQDIM],I,J));
