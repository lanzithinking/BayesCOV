% optimization for MCMC initialization
% for multiple trials

function [gamma,eta,mu,tau,vecL,objf]=opt4ini(gamma,eta,mu,tau,vecL,y,comm_ND,M_vecL,ker,a,b,m,V,objf,Nmax,thld)
addpath('../../optimizer/');
if ~exist('Nmax','var')
    Nmax=100;
end
if ~exist('thld','var')
    thld=1e-3;
end
if ~exist('objf','var') || isempty(objf)
    objf=inf(1,5);
end

% dimension
[M,N,sqdim]=size(y);

% optimization setting
h=1e-1; L=5;
r=.5; c=.1;

% define Dlt-SphOpt optimizer
dlt_sphOpt=Deltat_SphOpt(vecL,@(q,geord)geom_L(q,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),geord),h,L,r,c);

fprintf('Optimizing parameters...\n');
prog=0.05:0.05:1;
tic;
for iter=1:Nmax
    % record current value
    gamma_=gamma;
    eta_=eta;
    mu_=mu; % (N,sqdim)
    tau_=tau; % (N,sqdim)
    vecL_=vecL; % (N,dim)
    objf_=objf;
    
    % display the progress
    if ismember(iter,floor(Nmax.*prog))
        fprintf('%.0f%% of max iterations completed.\n',100*iter/Nmax);
    end
    
    if iter==0
        % update gamma
        % gamma_mu
        quad_mu=mu.*((exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(1)))+ker.jit)\mu);
        tr_mu=sum(quad_mu(:));
    %     gamma(1)=1./gamrnd(a(1)+N*sqdim/2,1/(b(1)+.5*tr_mu));
        alpha(1)=a(1)+N*sqdim/2; beta(1)=b(1)+.5*tr_mu;
        % gamma_tau
        quad_tau=tau.*((exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(2)))+ker.jit)\tau);
        tr_tau=sum(quad_tau(:));
    %     gamma(2)=1./gamrnd(a(2)+N*sqdim/2,1/(b(2)+.5*tr_tau));
        alpha(2)=a(2)+N*sqdim/2; beta(2)=b(2)+.5*tr_tau;
        % gamma_L
        quad_L=(vecL-M_vecL).*((exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(3)))+ker.jit)\(vecL-M_vecL));
        tr_L=sum(quad_L(:));
    %     gamma(3)=1./gamrnd(a(3)+N*sqdim/2*((sqdim+1)/2-1/sqdim),1/(b(3)+.5*tr_L));
        alpha(3)=a(3)+N*sqdim/2*((sqdim+1)/2-1/sqdim); beta(3)=b(3)+.5*tr_L;

        gamma=beta./(alpha+1);
        objf(1)=-sum(log(gampdf(1./gamma,alpha,1./beta)));

        % update eta
        % eta_mu
        logf_mu=@(q)logpost_eta(q,ker,gamma(1),mu,m(1),V(1));
    %     [eta(1),l_eta(1)]=slice(eta(1),logf_mu(eta(1)),logf_mu);
        % eta_tau
        logf_tau=@(q)logpost_eta(q,ker,gamma(2),tau,m(2),V(2));
    %     [eta(2),l_eta(2)]=slice(eta(2),logf_tau(eta(2)),logf_tau);
        % eta_L
        logf_L=@(q)logpost_eta(q,ker,gamma(3),vecL-M_vecL,m(3),V(3),sqdim*((sqdim+1)/2-1/sqdim));
    %     [eta(3),l_eta(3)]=slice(eta(3),logf_L(eta(3)),logf_L);
        [eta,objf(2)]=fminunc(@(q)-logf_mu(q(1))-logf_tau(q(2))-logf_L(q(3)),eta);
    end
    
    % update mu
    Rho_mat=zeros(sqdim,sqdim,N);
    for i=1:sqdim
        for j=1:sqdim
            ind_i=1+i*(i-1)/2:i*(i+1)/2;
            ind_j=1+j*(j-1)/2:j*(j+1)/2;
            min_ij=min([i,j]);
            Rho_mat(i,j,:)=sum(vecL(:,ind_i(1:min_ij)).*vecL(:,ind_j(1:min_ij)),2);
        end
    end
    sigma=exp(tau);
    Sigma_cell=num2cell(reshape(sigma',[sqdim,1,N]).*Rho_mat.*reshape(sigma',[1,sqdim,N]),[1,2]);
    Sigma_bkdg=sparse(blkdiag(Sigma_cell{:}));
    Sigma_perm=comm_ND'*Sigma_bkdg*comm_ND;
    Sigma_I=Sigma_perm*kron(speye(sqdim),inv(gamma(1).*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(1)))+ker.jit))+M.*speye(N*sqdim);
    Sigma_post=Sigma_I\Sigma_perm;
%     mu_vec=mvnrnd((Sigma_I\sum(y(:,:),1)')',(Sigma_post+Sigma_post')./2);
    mu_vec=(Sigma_I\sum(y(:,:),1)')';
    objf(3)=sum(log(diag(chol((Sigma_post+Sigma_post')./2))));
    mu=reshape(mu_vec,N,sqdim);
    
    % update tau
%     logLik_tau=@(q)loglik_tau(q,mu,y,Rho_mat,comm_ND);
    logf_tau=@(q)geom_tau(q,mu,y,Rho_mat,comm_ND,ker,gamma(2),eta(2));
%     rnd_pri_tau=@()mvnrnd(zeros(1,N),gamma(2).*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(2)))+ker.jit,sqdim)';
%     [tau,l_tau] = ESS(tau,logLik_tau(tau),rnd_pri_tau,logLik_tau);
%     K_tau=gamma(2).*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(2)))+ker.jit;
%     [tau,objf(4)]=fminunc(@(q)-logLik_tau(q)+sum(sum((K_tau\q).*q))./2,tau);
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'display','off',...
                           'MaxIterations',100);
%     options = optimset('Algorithm','trust-region','GradObj','on','Display','off');
    [tau,objf(4)]=fminunc(@(q)logf_tau(q),tau,options);
    
    % update L
    % update geometry function for L
    dlt_sphOpt.geom=@(q,geord)geom_L(q,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),geord);
    % update potential and its gradient for L
    [dlt_sphOpt.U,dU]=dlt_sphOpt.geom(vecL,[0,1]);
    dlt_sphOpt.DU=dlt_sphOpt.mult_proj(dU,vecL);
    dlt_sphOpt=dlt_sphOpt.backtrack(10);
%     dlt_sphOpt=dlt_sphOpt.optimize(10);
    vecL=dlt_sphOpt.q;
    objf(5)=dlt_sphOpt.U;
    
    % display current objective function
    fprintf(['Objective function values: ',repmat('%.4f, ',1,length(objf)),'at iteration %d.\n'], objf, iter);
    
    % break if condition satisfied
    if iter>1
        dif(1)=max(abs(gamma_-gamma)); dif(2)=max(abs(eta_-eta));
        dif(3)=max(abs(mu_(:)-mu(:))); dif(4)=max(abs(tau_(:)-tau(:)));
        dif(5)=max(abs(vecL_(:)-vecL(:)));
        if all(abs(objf_-objf)<thld) || all(dif<thld)
            fprintf('Optimization breaks at iteration %d.\n',iter);
            break;
        end
    end

end

% count time
time=toc;
fprintf('Time used: %.2f seconds.\n',time);
