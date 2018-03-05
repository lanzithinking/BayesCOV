% generate sample(s) from a unit multivariate Normal distribution

function r=umvnrnd(mu,Sigma,n)
[sz1,d]=size(mu);
if ~exist('Sigma','var')
    Sigma=eye(d);
end
if ~exist('n','var')
    n=sz1;
end
r=mvnrnd(mu(mod((1:n)-1,sz1)+1,:),Sigma,n);
r=r./sqrt(sum(r.^2,2));
end