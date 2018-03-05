% Multiple subscripts from linear index in vech
% ---------------------------------------------
% Shiwei Lan @ CalTech, 2017
% __________________________

function [I,J] = ind_vech2sub(SQDIM,IND,order)
if ~exist('order','var')
    order='col';
end

idx_mat = ivech(1:SQDIM*(SQDIM+1)/2,order);
ind_nz = idx_mat(:)~=0;
ref_idx = (1:SQDIM^2)';
ref_idx = ref_idx(ind_nz);
if any(strfind('row',order))
    ref_idx(idx_mat(ind_nz)) = ref_idx;
end
[I,J] = ind2sub(SQDIM,ref_idx(IND));




