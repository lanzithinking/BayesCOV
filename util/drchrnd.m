% generate sample(s) from a dirichlet distribution

function r = drchrnd(alpha,n)
sz1=size(alpha,1);
if ~exist('n','var')
    n=sz1;
end
r = gamrnd(alpha(mod((1:n)-1,sz1)+1,:),1);
r = r ./ sum(r,2);
end