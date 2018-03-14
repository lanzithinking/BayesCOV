% log-likelihood function of tau

function l=loglik_tau(tau,mu,y,Rho_mat,comm_ND)

isigma=exp(-tau);
y_std=(y-mu).*isigma;
Rho_cell=num2cell(Rho_mat,[1,2]);
Rho_bkdg=sparse(blkdiag(Rho_cell{:}));
Rho_perm=comm_ND'*Rho_bkdg*comm_ND;
l=-sum(tau(:))-.5*y_std(:)'/Rho_perm*y_std(:);

end