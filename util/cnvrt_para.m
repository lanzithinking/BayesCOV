% convert samples of tau, L to sigma, Rho

function [sigma,vecRho]=cnvrt_para(tau,vecL,rev)
if ~exist('rev','var')
    rev=false;
end

sigma=exp(tau);

L=ivech(vecL,'row');
if rev
    d=(sqrt(8*length(vecL)+1)-1)/2;
    L=L(d:-1:1,d:-1:1); % convert it to upper triangular matrix
end
vecRho=vech(L*L','row');
