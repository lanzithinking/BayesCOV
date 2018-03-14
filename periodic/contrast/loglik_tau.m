% log-likelihood function of tau
% for multiple trials

function l=loglik_tau(tau,mu,y,Rho_mat,comm_ND)

isigma=exp(-tau);
y_std=(y-shiftdim(mu,-1)).*shiftdim(isigma,-1); % (M,N,sqdim)
Rho_cell=num2cell(Rho_mat,[1,2]);
Rho_bkdg=sparse(blkdiag(Rho_cell{:}));
Rho_perm=comm_ND'*Rho_bkdg*comm_ND;
yinvR=y_std(:,:)/Rho_perm; % (M,Nsqdim)
l=-size(y,1)*sum(tau(:))-.5*sum(sum(yinvR.*y_std(:,:)));

end