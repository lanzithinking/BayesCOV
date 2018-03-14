% inference of a periodic process model using Dlt-Opt and Dlt-SphHMC
% for multiple-trials

clear;
addpath('../','../../util/','../../sampler/','../../optimizer/');
% Random Numbers...
seedNO=2017;
seed=RandStream('mt19937ar','Seed',seedNO);
RandStream.setGlobalStream(seed);

% setting for simulation
M=20; N=100; sqdim=10;
% load or simulate data
[t,y]=generate_data(M,N,sqdim,seedNO);
dim=sqdim*(sqdim+1)/2;
% constant matrix
comm_ND=sparse(vecperm(N,sqdim)');

% mean in uvGP for vecL
M_=eye(sqdim);
% % M_=chol(corr(reshape(y,[],sqdim)),'lower');
% M_=corr(reshape(y,[],sqdim));
% M_=sign(M_).*double(abs(M_)>0.9);M_(1:sqdim+1:end)=1;
M_vecL=vech(M_,'row')'; % (1,dim)
% parameters of kernel
s=2; % smoothness
dist_t=pdist2(t,t,@(XI,XJ)abs(bsxfun(@minus,XI,XJ)).^s);
jit=1e-5.*speye(N);
ker.s=s;ker.dist_t=dist_t;ker.jit=jit;
% specify (hyper)-priors
a=ones(1,3); b=[.1,1e-3,2.5e-1]; % (a,b) in inv-gamma priors for gamma_*, * = mu, tau(log-sigma), L
m=zeros(1,3); V=[1,1,2]; % (m,V) in (log) normal priors for eta_*, (eta=log-rho), * = mu, tau(log-sigma), L
% all: a=ones(1,3); m=zeros(1,3);
% N=20: b=[.1,1e-3,.2]; V=[1,.5,1];
% N=200: b=[.1,1e-3,.2]; V=[1,1,.3];

% optimize gamma, eta
% sample mu, tau, L using Gibbs, ESS, Dlt-SphHMC, resp.

% sampling setting
stepsz=1e-3; Nleap=20;
alg={'standard','robust'};
alg_choice=alg{1};
reweight=false;

% allocation to save
Nsamp=1e3; burnrate=1; thin=10;
Niter=Nsamp*thin; NBurnIn=floor(Niter*burnrate); Niter=Niter+ NBurnIn;
samp_mu=zeros(Nsamp,N,sqdim);
samp_tau=zeros(Nsamp,N,sqdim);
samp_vecL=zeros(Nsamp,N,dim);
samp_gamma=zeros(Nsamp,3);
samp_eta=zeros(Nsamp,3);
engy=zeros(Niter,5);
l_eta=zeros(1,3);
acpt=0; % overall acceptance
accp=0; % online acceptance

% initializatioin
gamma=1./gamrnd(a,1./b);
eta=normrnd(m,sqrt(V));
mu=mvnrnd(zeros(1,N),gamma(1).*exp(-.5.*dist_t.*exp(-s.*eta(1)))+jit,sqdim)'; % (N,sqdim)
tau=mvnrnd(zeros(1,N),gamma(2).*exp(-.5.*dist_t.*exp(-s.*eta(2)))+jit,sqdim)'; % (N,sqdim)
vecL=mvnrnd(M_vecL'*ones(1,N),gamma(3).*exp(-.5.*dist_t.*exp(-s.*eta(3)))+jit)'; % (N,dim)
for d=1:sqdim
    idx_d=1+d*(d-1)/2:d*(d+1)/2;
    vecL(:,idx_d)=vecL(:,idx_d)./sqrt(sum(vecL(:,idx_d).^2,2));
end
vecL(:,1)=1;

options=optimoptions('fminunc','Algorithm','quasi-newton','display','off','MaxIterations',100);
% optimize initial location
[gamma,eta,mu,tau,vecL,objf]=opt4ini(gamma,eta,mu,tau,vecL,y,comm_ND,M_vecL,ker,a,b,m,V,[],20);

% define Dlt-SphHMC sampler
dlt_sphHMC4L=Deltat_SphHMC(vecL,@(q,geord,adj_metric)geom_L(q,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),geord,adj_metric),stepsz,Nleap,alg_choice,reweight);

fprintf('Running %s Delta-Spherical HMC with %s reweight...\n',alg_choice,~reweight*'no');
prog=0.05:0.05:1;
tic;
for iter=1:Niter
    
    % display sampleing progress and online acceptance rate
    if ismember(iter,floor(Niter.*prog))
        fprintf('%.0f%% iterations completed.\n',100*iter/Niter);
        fprintf('Online acceptance rate: %.0f%%\n', accp.*(100/floor(prog(1)*Niter)));
        accp=0;
    end
    
    % optimize for gamma
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
    % update gamma
    gamma=beta./(alpha+1);
    
    % optimize for eta
    % eta_mu
    logf_mu=@(q)logpost_eta(q,ker,gamma(1),mu,m(1),V(1));
%     [eta(1),l_eta(1)]=slice(eta(1),logf_mu(eta(1)),logf_mu);
    [eta(1),l_eta(1)]=fminunc(@(q)-logf_mu(q),eta(1),options);
    % eta_tau
    logf_tau=@(q)logpost_eta(q,ker,gamma(2),tau,m(2),V(2));
%     [eta(2),l_eta(2)]=slice(eta(2),logf_tau(eta(2)),logf_tau
    [eta(2),l_eta(2)]=fminunc(@(q)-logf_tau(q),eta(2),options);
    % eta_L
    logf_L=@(q)logpost_eta(q,ker,gamma(3),vecL-M_vecL,m(3),V(3),sqdim*((sqdim+1)/2-1/sqdim));
%     [eta(3),l_eta(3)]=slice(eta(3),logf_L(eta(3)),logf_L);
    [eta(3),l_eta(3)]=fminunc(@(q)-logf_L(q),eta(3),options);
    % update eta
%     [eta,l_eta]=fminunc(@(q)-logf_mu(q(1))-logf_tau(q(2))-logf_L(q(3)),eta,options);
    
    % sample mu
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
    mu_vec=mvnrnd((Sigma_I\sum(y(:,:),1)')',(Sigma_post+Sigma_post')./2);
    mu=reshape(mu_vec,N,sqdim);
    
    % sample tau with ESS
    logLik_tau=@(q)loglik_tau(q,mu,y,Rho_mat,comm_ND);
    rnd_pri_tau=@()mvnrnd(zeros(1,N),gamma(2).*exp(-.5.*ker.dist_t.*exp(-ker.s.*eta(2)))+ker.jit,sqdim)';
    [tau,l_tau] = ESS(tau,logLik_tau(tau),rnd_pri_tau,logLik_tau);
    
    % sample L with Delta-Spherical HMC
    % update geometry function for L
    dlt_sphHMC4L.geom=@(q,geord)geom_L(q,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),geord,~dlt_sphHMC4L.reweight);
    % update potential and its gradient for L
    [dlt_sphHMC4L.U,dlt_sphHMC4L.dU]=dlt_sphHMC4L.geom(vecL,[0,1]);
    switch alg_choice
        case alg{1}
            [dlt_sphHMC4L,acpt_L]=dlt_sphHMC4L.std_Dlt_SphHMC(true);
        case alg{2}
            [dlt_sphHMC4L,acpt_L]=dlt_sphHMC4L.rbt_Dlt_SphHMC(true);
    end
    vecL=dlt_sphHMC4L.q;
    
    % accptances
    accp=accp+acpt_L;
    
    % burn-in complete
    if(iter==NBurnIn)
        fprintf('Burn in completed!\n');
    end
    
    engy(iter,:)=[l_eta,l_tau,dlt_sphHMC4L.U];
    % save samples after burn-in
    if(iter>NBurnIn)
        if mod(iter-NBurnIn,thin) == 0
            samp_gamma(ceil((iter-NBurnIn)/thin),:)=gamma;
            samp_eta(ceil((iter-NBurnIn)/thin),:)=eta;
            samp_mu(ceil((iter-NBurnIn)/thin),:,:)=mu;
            samp_tau(ceil((iter-NBurnIn)/thin),:,:)=tau;
            samp_vecL(ceil((iter-NBurnIn)/thin),:,:)=vecL;
        end
        acpt=acpt+acpt_L;
    end

end

% count time
time=toc;
% save results
acpt=acpt/(Niter-NBurnIn);
% resample if set to reweight
if reweight
    samp_diag=samp_vecL(:,:,(1:sqdim).*((1:sqdim)+1)./2);
    wt=sum(log(abs(samp_diag(:,:))),2);
    wt=exp(wt-median(wt));
    resamp_idx=datasample(1:Nsamp,Nsamp,'replace',true,'weights',wt);
else
    wt=ones(Nsamp,1);
    resamp_idx=1:Nsamp;
end
% save
time_lbl=regexprep(num2str(fix(clock)),'    ','_');
if reweight
    ifrewt='_reweight';
else
    ifrewt='_noreweight';
end
f_name=['periodic_dlt_OptSphHMC_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim),'_',time_lbl,ifrewt];
if exist('result','dir')~=7
    mkdir('result');
end
save(['result/',f_name,'.mat'],'seedNO','t','y','M_vecL','ker','a','b','m','V',...
                               'alg_choice','Nsamp','NBurnIn','thin','stepsz','Nleap',...
                               'samp_gamma','samp_eta','samp_mu','samp_tau','samp_vecL','reweight','wt','resamp_idx','acpt','time');
% summarize
fprintf('\nIt takes %.2f seconds to collect %e samples after thinning %d.\n', time,Nsamp,thin);
fprintf('\nThe final acceptance rate is: %.4f.\n',acpt);
% efficiency measurement
addpath('./result/');
% CalculateStatistics(f_name,'./result/');