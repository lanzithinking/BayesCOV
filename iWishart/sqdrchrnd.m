 % generate sample(s) from a squared dirichlet distribution

function r = sqdrchrnd(a,n)
if ~exist('n','var')
    n=1;
end
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);
s = 1 - 2.*binornd(1,0.5,n,p);
r = sqrt(r) .* s;
end