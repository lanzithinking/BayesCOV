% log-posterior density of eta

function l=logpost_eta(eta,ker,gamma,star,m,V,sz)
if ~exist('sz','var')
    sz=size(star,2);
end

chol_K=chol(gamma.*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta))+ker.jit,'lower');
half_quad=chol_K\star;
l=-sz*sum(log(diag(chol_K)))-.5*sum(half_quad(:).^2)-.5*(eta-m).^2./V;

end