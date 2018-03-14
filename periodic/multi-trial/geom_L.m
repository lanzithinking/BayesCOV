% Energy function for cholesky factor L %
% for multiple trials %

function [U,dU]=geom_L(vecL,mu,tau,y,comm_ND,M_vecL,ker,gamma,eta,geord,adj_metric)
U=[];dU=[];
if ~exist('geord','var')
    geord=0;
end
if ~exist('adj_metric','var')
    adj_metric=false;
end

[M,N,sqdim]=size(y); dim=size(vecL,2);
% for likelihood
isigma=exp(-tau);
y_std=(y-shiftdim(mu,-1)).*shiftdim(isigma,-1); % (M,N,sqdim)
% for prior
K_L=gamma.*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta))+ker.jit;
invKl=K_L\(vecL-M_vecL); % (N,dim)

if any(geord==0)
%     tic;
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
    yinvR=y_std(:,:)/Rho_perm; % (M,Nsqdim)
    loglik=-sum(log(abs(diagL(:)))).*(M+adj_metric)-.5*sum(sum(yinvR.*y_std(:,:)));
%     time=toc;
%     fprintf('Time used: %.4f\n',time);
%     tic;
%     L_mat=ivechx(vecL,'row');
%     L_cell=num2cell(permute(L_mat,[2,3,1]),[1,2]);
%     L_bkdg=sparse(blkdiag(L_cell{:}));
%     L_perm=comm_ND'*L_bkdg;
%     yinvL=y_std(:,:)/L_perm'; % (M,Nsqdim)
%     loglik1=-sum(log(abs(diagL(:)))).*(M+adj_metric)-.5*sum(sum(yinvL.^2));
%     time=toc;
%     fprintf('Time used: %.4f\n',time);
%     fprintf('Difference: %.8f\n',loglik-loglik1);

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
        yinvL_n=reshape((y_std(:,n,:)),M,sqdim)/L_n'; % (M,sqdim)
        dloglik(n,:)=vech(-diag(1./diag(L_n)).*(M+adj_metric)+L_n'\yinvL_n'*yinvL_n,'row');
    end
%     dloglik1=spdiags(-(M+adj_metric).*spdiags(L_bkdg),0,size(L_bkdg)) + L_bkdg'\yinvL'*yinvL; % then extract block diagonals: not worth it
    
    % gradient of log-prior
    dlogpri=-invKl; % (N,dim)
    dU=-(dloglik+dlogpri);
end

end