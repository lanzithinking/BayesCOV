% Energy function of tau

function [U,dU]=geom_tau(tau,mu,y,Rho_mat,comm_ND,ker,gamma,eta,geord)
U=[];dU=[];
if ~exist('geord','var')
    geord=0:1;
end

% for likelihood
isigma=exp(-tau);
y_std=(y-mu).*isigma;
Rho_cell=num2cell(Rho_mat,[1,2]);
Rho_bkdg=sparse(blkdiag(Rho_cell{:}));
Rho_perm=comm_ND'*Rho_bkdg*comm_ND;
% loglik=-sum(tau(:))-.5*y_std(:)'/Rho_perm*y_std(:);
invRhoy=Rho_perm\y_std(:);
% for prior
K_tau=gamma.*(exp(-.5.*ker.dist_t.*exp(-ker.s.*eta))+ker.jit);
invKtau=K_tau\tau; % (N,sqdim)

if any(geord==0)
    % log-likelihood
    loglik=-sum(tau(:))-.5*y_std(:)'*invRhoy;
    % log-prior
    quad_tau=tau.*invKtau;
    logpri=-.5.*sum(quad_tau(:));
    
    U=-(loglik+logpri);
end

if any(geord==1)
    % gradient of log-likelihood
    dloglik=reshape(-1+y_std(:).*invRhoy,size(y));
    % gradient of log-prior
    dlogpri=-invKtau; % (N,sqdim)
    
    dU=-(dloglik+dlogpri);
end