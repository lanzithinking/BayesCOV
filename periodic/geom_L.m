% Energy function for cholesky factor L %

function [U,dU]=geom_L(vecL,mu,tau,y,comm_ND,M_vecL,ker,gamma,eta,geord,adj_metric)
U=[];dU=[];
if ~exist('geord','var')
    geord=0;
end
if ~exist('adj_metric','var')
    adj_metric=false;
end

[N,dim]=size(vecL); sqdim=size(mu,2);
% for likelihood
isigma=exp(-tau);
y_std=(y-mu).*isigma;
% for prior
K_L=gamma.*(exp(-.5.*ker.dist_t.*exp(-ker.s.*eta))+ker.jit);
invKl=K_L\(vecL-M_vecL); % (N,dim)

if any(geord==0)
    Rho_mat=zeros(sqdim,sqdim,N);
    for i=1:sqdim
        for j=1:sqdim
            ind_i=1+i*(i-1)/2:i*(i+1)/2;
            ind_j=1+j*(j-1)/2:j*(j+1)/2;
            min_ij=min([i,j]);
            Rho_mat(i,j,:)=sum(vecL(:,ind_i(1:min_ij)).*vecL(:,ind_j(1:min_ij)),2);
        end
    end
    Rho_cell=num2cell(Rho_mat,[1,2]);
    Rho_bkdg=sparse(blkdiag(Rho_cell{:}));
    Rho_perm=comm_ND'*Rho_bkdg*comm_ND;
    % log-likelihood
    diagL=vecL(:,(1:sqdim).*((1:sqdim)+1)./2);
    loglik=-sum(log(abs(diagL(:)))).*(1+adj_metric)-.5*y_std(:)'/Rho_perm*y_std(:);

    quad_L=(vecL-M_vecL).*invKl;
    % log-prior
    logpri=-.5.*sum(quad_L(:));

    U=-(loglik+logpri);
end

if any(geord==1)
    % gradient of log-likelihood
    dloglik=zeros(N,dim);
    for n=1:N
        L_n=ivech(vecL(n,:),'row');
        invLy_n=L_n\y_std(n,:)';
%         dloglik(n,:)=vech(-diag(1./diag(L_n)).*(1+adj_metric)+tril(L_n'\invLy_n*invLy_n'),'row');
        dloglik(n,:)=vech(-diag(1./diag(L_n)).*(1+adj_metric)+L_n'\invLy_n*invLy_n','row');
    end
    % gradient of log-prior
    dlogpri=-invKl; % (N,dim)
    dU=-(dloglik+dlogpri);
end

end