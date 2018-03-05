% log-likelihood function %

function l=loglik(tau,vecL,prior,data)

isigma=exp(-tau);
y_std=(data.y'-data.mu0).*isigma;

L=ivech(vecL,'row');
if contains(prior.choice,'bartlett')
    L=L(end:-1:1,end:-1:1); % becomes upper triangular matrix
end
diagL=diag(L);
invLy_std=L\y_std;
l=-data.size*(sum(tau)+sum(log(abs(diagL))))-sum(invLy_std(:).^2)./2;

end