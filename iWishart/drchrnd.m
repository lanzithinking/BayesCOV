% generate sample(s) from a dirichlet distribution

function r = drchrnd(a,n)
if ~exist('n','var')
    n=1;
end
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);
end