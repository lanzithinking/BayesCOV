 % generate sample(s) from a squared dirichlet distribution

function r = sqdrchrnd(alpha,n)
[sz1,d]=size(alpha);
if ~exist('n','var')
    n=sz1;
end
r = gamrnd(alpha(mod((1:n)-1,sz1)+1,:),1);
r = r ./ sum(r,2);
s = 1 - 2.*binornd(1,0.5,n,d);
r = sqrt(r) .* s;
end