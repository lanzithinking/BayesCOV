% Energy function of tau
% for multiple trials

function [U,dU]=geom_tau(tau,mu,y,Rho_mat,comm_ND,ker,gamma,eta,geord)
U=[];dU=[];
if ~exist('geord','var')
    geord=0:1;
end

M=size(y,1);
% for likelihood
isigma=exp(-tau);
y_std=(y-shiftdim(mu,-1)).*shiftdim(isigma,-1); % (M,N,sqdim)
Rho_cell=num2cell(Rho_mat,[1,2]);
Rho_bkdg=sparse(blkdiag(Rho_cell{:}));
Rho_perm=comm_ND'*Rho_bkdg*comm_ND;
% yinvR=y_std(:,:)/Rho_perm; % (M,Nsqdim)
% quad_y=sum(yinvR.*y_std(:,:),1); % (1,Nsqdim)
quad_y=sum((y_std(:,:)/Rho_perm).*y_std(:,:),1); % (1,Nsqdim)
% for prior
K_tau=gamma.*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta))+ker.jit;
invKtau=K_tau\tau; % (N,sqdim)

if any(geord==0)
    % log-likelihood
    loglik=-M*sum(tau(:))-.5*sum(quad_y);
    % log-prior
    quad_tau=tau.*invKtau;
    logpri=-.5.*sum(quad_tau(:));
    
    U=-(loglik+logpri);
end

if any(geord==1)
    % gradient of log-likelihood
    dloglik=reshape(-M+quad_y,size(tau));
    % gradient of log-prior
    dlogpri=-invKtau; % (N,sqdim)
    
    dU=-(dloglik+dlogpri);
end